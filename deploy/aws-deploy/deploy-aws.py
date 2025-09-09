#!/usr/bin/env python3
"""
AWS Pool-Based Deployment Script for BioML-bench

Manages a persistent pool of EC2 instances to run biomlbench jobs:
1. Creates/uses a pool of EC2 instances with biomlbench AMI
2. Pre-warms instances for optimal performance
3. Runs biomlbench pipeline using dynamic job queue
4. Verifies S3 uploads
5. Keeps instances running for reuse (manual cleanup required)
"""

import argparse
import subprocess
import time
import uuid
import json
import boto3
import queue
import threading
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple, Optional

# AWS clients
ec2 = boto3.client('ec2')

# Default configuration
DEFAULT_REGION = "us-east-1"
DEFAULT_KEY_NAME = "biomlbench-key"  # Matches setup-aws-resources.sh

def log(message: str):
    """Simple logging with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def run_command(cmd: List[str], description: str = "", timeout: int = None) -> Tuple[int, str]:
    """Run a command and return (exit_code, output)."""
    try:
        if description:
            log(f"Running: {description}")
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            timeout=timeout
        )
        return result.returncode, result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return 124, "Command timed out"
    except Exception as e:
        return 1, str(e)

def get_instance_type(machine_type: str) -> str:
    """Map machine type choice to AWS instance type."""
    if machine_type == "cpu":
        # n2-standard-16 (GCP) → m5.4xlarge (AWS)
        # 16 vCPUs, 64GB RAM → 16 vCPUs, 64GB RAM
        return "m5.4xlarge"
    else:  # gpu
        # g2-standard-16 (GCP with L4) → g4dn.4xlarge (AWS with T4)
        # Similar GPU performance tier
        return "g4dn.4xlarge"

def create_instance_with_retry(instance_name: str, machine_type: str, ami_id: str, 
                              key_name: str = DEFAULT_KEY_NAME,
                              security_group: str = DEFAULT_SECURITY_GROUP,
                              max_retries: int = None,
                              custom_tag: str = None) -> Optional[str]:
    """Create an EC2 instance, retrying until success. Returns instance ID if successful."""
    log(f"Creating instance: {instance_name} (machine type: {machine_type})")
    
    instance_type = get_instance_type(machine_type)
    
    # User data script to set up the instance
    user_data = """#!/bin/bash
# Start Docker daemon if not running
systemctl start docker || true

# Ensure runner user can access docker
usermod -aG docker runner || true
chmod 666 /var/run/docker.sock || true
"""
    
    attempt = 0
    while True:
        attempt += 1
        if max_retries and attempt > max_retries:
            log(f"❌ Failed to create {instance_name} after {max_retries} attempts")
            return None
        
        try:
            response = ec2.run_instances(
                ImageId=ami_id,
                InstanceType=instance_type,
                MinCount=1,
                MaxCount=1,
                KeyName=key_name,
                SecurityGroupIds=[security_group],
                BlockDeviceMappings=[
                    {
                        'DeviceName': '/dev/sda1',  # or /dev/xvda depending on AMI
                        'Ebs': {
                            'VolumeSize': 1000,
                            'VolumeType': 'gp3',
                            'DeleteOnTermination': True
                        }
                    }
                ],
                UserData=user_data,
                TagSpecifications=[
                    {
                        'ResourceType': 'instance',
                        'Tags': [
                            {'Key': 'Name', 'Value': instance_name},
                            {'Key': 'Project', 'Value': 'biomlbench'},
                            {'Key': 'ManagedBy', 'Value': 'deploy-aws.py'}
                        ] + ([{'Key': 'CustomTag', 'Value': custom_tag}] if custom_tag else [])
                    }
                ]
            )
            
            instance_id = response['Instances'][0]['InstanceId']
            log(f"✅ Created instance: {instance_name} (ID: {instance_id})")
            
            # Wait for instance to be running
            ec2.get_waiter('instance_running').wait(InstanceIds=[instance_id])
            
            return instance_id
            
        except Exception as e:
            log(f"⚠️  Instance creation failed (attempt {attempt}): {str(e)}")
            time.sleep(5)  # Wait before retry

def get_instance_ip(instance_id: str) -> Optional[str]:
    """Get the public IP of an instance."""
    try:
        response = ec2.describe_instances(InstanceIds=[instance_id])
        public_ip = response['Reservations'][0]['Instances'][0].get('PublicIpAddress')
        return public_ip
    except Exception as e:
        log(f"❌ Failed to get IP for {instance_id}: {str(e)}")
        return None

def wait_for_ssh(instance_id: str, timeout: int = 300) -> bool:
    """Wait for SSH to become available on the instance."""
    log(f"Waiting for SSH on instance {instance_id}...")
    
    public_ip = get_instance_ip(instance_id)
    if not public_ip:
        return False
    
    log(f"Instance IP: {public_ip}")
    
    start_time = time.time()
    while time.time() - start_time < timeout:
        exit_code, output = run_command([
            "ssh",
            "-o", "StrictHostKeyChecking=no",
            "-o", "ConnectTimeout=10",
            "-o", "BatchMode=yes",
            f"runner@{public_ip}",
            "echo 'SSH ready'"
        ], timeout=15)
        
        if exit_code == 0:
            log(f"✅ SSH ready on {public_ip}")
            return True
        else:
            log(f"SSH attempt failed: {output.strip()}")
        
        time.sleep(10)
    
    log(f"❌ SSH timeout on instance {instance_id}")
    return False

def run_biomlbench_job(instance_id: str, agent: str, task_id: str) -> bool:
    """Run the biomlbench pipeline on an instance using SSH."""
    log(f"Running job on {instance_id}: {agent} -> {task_id}")
    
    public_ip = get_instance_ip(instance_id)
    if not public_ip:
        return False
    
    # Build the command to run on the instance
    remote_cmd = f"""
    set -e
    cd /home/runner/biomlbench
    source .venv/bin/activate
    
    # Run agent with fast container config
    biomlbench run-agent --s3-bucket biomlbench --s3-prefix v3/artifacts --agent {agent} --task-id {task_id} --cpu-only --container-config environment/config/container_configs/fast.json
    
    # Get the specific run group ID that was just created
    LATEST_RUN=$(find runs/ -name "*run-group_{agent}" -type d | sort | tail -1)
    RUN_GROUP_ID=$(basename "$LATEST_RUN")
    echo "📁 Run group: $RUN_GROUP_ID"
    
    # Grade results
    biomlbench grade --s3-bucket biomlbench --s3-prefix v3/artifacts --submission "$LATEST_RUN/submission.jsonl" --output-dir results/
    
    # Get the grading timestamp from the most recent grading report
    GRADING_REPORT=$(find results/ -name "*_grading_report.json" | sort | tail -1)
    GRADING_TIMESTAMP=$(basename "$GRADING_REPORT" | cut -d'_' -f1)
    echo "📊 Grading timestamp: $GRADING_TIMESTAMP"
    
    # Convert task_id to S3-safe format for the organized structure
    TASK_ID_SAFE=$(echo "{task_id}" | sed 's/\\//-/g' | sed 's/_/-/g')
    
    # Show the exact S3 paths for this specific run
    echo "📤 S3 artifacts for this run:"
    echo "  Run artifacts:"
    echo "    s3://biomlbench/v3/artifacts/runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz"
    echo "    OR s3://biomlbench/v3/artifacts/failed_runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz (if failed)"
    echo "  Grading artifacts:"
    echo "    s3://biomlbench/v3/artifacts/grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_grading_report.json.gz"
    echo "    s3://biomlbench/v3/artifacts/grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_individual_reports.tar.gz"
    echo "    OR s3://biomlbench/v3/artifacts/failed_grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_*.gz (if failed)"
    
    # Verify these specific paths exist
    echo "🔍 Verifying uploads..."
    if aws s3 ls s3://biomlbench/v3/artifacts/runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Run artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v3/artifacts/failed_runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Failed run artifacts uploaded successfully (organized structure)"
    else
        echo "❌ No run artifacts found in S3!"
        exit 1
    fi
    
    if aws s3 ls s3://biomlbench/v3/artifacts/grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_grading_report.json.gz > /dev/null 2>&1; then
        echo "✅ Grading artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v3/artifacts/grades/${{GRADING_TIMESTAMP}}_grading_report.json.gz > /dev/null 2>&1; then
        echo "✅ Grading artifacts uploaded successfully (flat structure)"
    elif aws s3 ls s3://biomlbench/v3/artifacts/failed_grades/{agent}/$TASK_ID_SAFE/ | grep -q "$GRADING_TIMESTAMP" > /dev/null 2>&1; then
        echo "✅ Failed grading artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v3/artifacts/failed_grades/ | grep -q "$GRADING_TIMESTAMP" > /dev/null 2>&1; then
        echo "✅ Failed grading artifacts uploaded successfully (flat structure)"
    else
        echo "❌ No grading artifacts found in S3!"
        exit 1
    fi
    
    echo "✅ Job completed successfully"
    """
    
    exit_code, output = run_command([
        "ssh",
        "-o", "StrictHostKeyChecking=no",
        f"runner@{public_ip}",
        remote_cmd
    ])
    
    if exit_code == 0:
        log(f"✅ Job completed on {instance_id}: {agent} -> {task_id}")
        # Show the output
        if output.strip():
            print("\n" + "="*60)
            print("JOB OUTPUT:")
            print("="*60)
            print(output.strip())
            print("="*60 + "\n")
        return True
    else:
        log(f"❌ Job failed on {instance_id}: {agent} -> {task_id}")
        log(f"Error output: {output}")
        return False





def load_jobs(jobs_file: str) -> List[Tuple[str, str]]:
    """Load jobs from a text file. Handles paths relative to current working directory."""
    jobs = []
    
    # Handle both absolute and relative paths
    jobs_path = Path(jobs_file)
    if not jobs_path.is_absolute():
        # If relative, try current working directory first, then relative to this script
        if (Path.cwd() / jobs_path).exists():
            jobs_path = Path.cwd() / jobs_path
        else:
            jobs_path = Path(__file__).parent / jobs_path
    
    try:
        with open(jobs_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                parts = line.split(',')
                if len(parts) != 2:
                    log(f"⚠️  Skipping invalid job on line {line_num}: {line}")
                    continue
                    
                agent, task_id = parts[0].strip(), parts[1].strip()
                jobs.append((agent, task_id))
        
        log(f"Loaded {len(jobs)} jobs from {jobs_path}")
        return jobs
    except FileNotFoundError:
        log(f"❌ Jobs file not found: {jobs_file}")
        log(f"   Looked for: {jobs_path}")
        return []

def verify_aws_setup(ami_id: str, key_name: str, security_group: str) -> bool:
    """Verify that all AWS resources exist."""
    try:
        # Check AMI
        response = ec2.describe_images(ImageIds=[ami_id])
        if not response['Images']:
            log(f"❌ AMI not found: {ami_id}")
            return False
        
        # Check key pair
        try:
            ec2.describe_key_pairs(KeyNames=[key_name])
        except ec2.exceptions.ClientError:
            log(f"❌ Key pair not found: {key_name}")
            return False
        
        # Check security group
        try:
            response = ec2.describe_security_groups(GroupIds=[security_group])
        except ec2.exceptions.ClientError:
            log(f"❌ Security group not found: {security_group}")
            return False
        
        log("✅ All AWS resources verified")
        return True
        
    except Exception as e:
        log(f"❌ Error verifying AWS setup: {str(e)}")
        return False

def create_instance_pool(count: int, machine_type: str, ami_id: str,
                        key_name: str, security_group: str, custom_tag: str = None) -> List[str]:
    """Create a pool of instances and return their IDs."""
    log(f"Creating pool of {count} instances...")
    instance_ids = []
    
    with ThreadPoolExecutor(max_workers=count) as executor:
        futures = []
        for i in range(count):
            instance_name = f"bioml-pool-{i+1}-{uuid.uuid4().hex[:8]}"
            future = executor.submit(
                create_instance_with_retry,
                instance_name, machine_type, ami_id,
                key_name, security_group, max_retries=3, custom_tag=custom_tag
            )
            futures.append(future)
        
        for future in as_completed(futures):
            instance_id = future.result()
            if instance_id:
                instance_ids.append(instance_id)
    
    log(f"Created {len(instance_ids)} instances successfully")
    return instance_ids

def prewarm_instances(instance_ids: List[str]):
    """Pre-warm critical paths on all instances using smart batching."""
    log(f"Pre-warming {len(instance_ids)} instances...")
    
    def warm_instance(instance_id):
        public_ip = get_instance_ip(instance_id)
        if not public_ip:
            return False
        
        # Wait for SSH first
        if not wait_for_ssh(instance_id):
            return False
        
        log(f"Pre-warming {instance_id} ({public_ip})...")
        
        # First, update biomlbench codebase to get latest changes
        log(f"Updating biomlbench codebase on {instance_id}...")
        
        update_cmd = """
        cd /home/runner/biomlbench
        git stash push -m "Stashing local changes before update" || true
        # Remove any conflicting files that might prevent git pull
        rm -f environment/config/container_configs/fast.json || true
        git pull origin main || git pull origin master || echo "Git pull failed, continuing with existing code"
        echo "Updated to commit: $(git rev-parse --short HEAD)"
        
        # Rebuild AIDE image
        echo "Rebuilding AIDE agent..."
        source .venv/bin/activate && bash scripts/build_agent.sh aide
        
        # Copy task descriptions to cache locations
        echo "Copying task descriptions to cache..."
        mkdir -p /home/runner/.cache/bioml-bench/data/manual/
        for task_dir in /home/runner/biomlbench/biomlbench/tasks/manual/*/; do
            if [ -d "$task_dir" ]; then
                task_name=$(basename "$task_dir")
                mkdir -p "/home/runner/.cache/bioml-bench/data/manual/$task_name/prepared/public/"
                if [ -f "$task_dir/description.md" ]; then
                    cp "$task_dir/description.md" "/home/runner/.cache/bioml-bench/data/manual/$task_name/prepared/public/description.md"
                    echo "Copied description for $task_name"
                else
                    echo "ERROR: Missing description.md for task $task_name"
                    exit 1
                fi
            fi
        done
        
        # Replace MLAgentBench LLM.py with custom version
        echo "Updating MLAgentBench LLM.py..."
        if [ -f "/home/runner/biomlbench/scripts/aws-deploy/LLM.py" ]; then
            cp "/home/runner/biomlbench/scripts/aws-deploy/LLM.py" "/home/runner/biomlbench/agents/mlagentbench/ref/MLAgentBench/MLAgentBench/LLM.py"
            echo "Replaced LLM.py with custom version"
        else
            echo "ERROR: Custom LLM.py not found at /home/runner/biomlbench/scripts/aws-deploy/LLM.py"
            exit 1
        fi
        
        # Rebuild MLAgentBench
        echo "Rebuilding MLAgentBench agent..."
        source .venv/bin/activate && bash scripts/build_agent.sh mlagentbench
        
        # Replace STELLA stella_core.py with custom version
        echo "Updating STELLA stella_core.py..."
        if [ -f "/home/runner/biomlbench/scripts/aws-deploy/stella_core.py" ]; then
            cp "/home/runner/biomlbench/scripts/aws-deploy/stella_core.py" "/home/runner/biomlbench/agents/stella/ref/STELLA/stella_core.py"
            echo "Replaced stella_core.py with custom version"
        else
            echo "ERROR: Custom stella_core.py not found at /home/runner/biomlbench/scripts/aws-deploy/stella_core.py"
            exit 1
        fi
        
        # Rebuild STELLA
        echo "Rebuilding STELLA agent..."
        source .venv/bin/activate && bash scripts/build_agent.sh stella
        
        echo "All updates completed successfully"
        """
        
        ssh_update_cmd = [
            "ssh", "-o", "StrictHostKeyChecking=no",
            f"runner@{public_ip}",
            update_cmd
        ]
        subprocess.run(ssh_update_cmd, check=False)  # Don't fail if git pull fails
        
        # Create fast entrypoint and config files on remote instance
        log(f"Creating fast configuration files on {instance_id}...")
        
        setup_cmd = """
        # Create fast entrypoint script
        cat > /home/runner/fast-entrypoint.sh << 'EOF'
#!/bin/bash

# Print commands and their arguments as they are executed
set -x

{
  # log into /home/logs
  LOGS_DIR=/home/logs
  mkdir -p $LOGS_DIR

  echo "Starting FAST entrypoint.sh at $(date)"
  echo "Directory contents before setup:"
  ls -la /home/

  # FAST permission setup - only what's absolutely necessary
  echo "Starting minimal permission setup at $(date)"
  
  # Make /home accessible
  chmod 755 /home 2>/dev/null || true
  
  # Only set up directories that actually exist and are needed
  if [ -d "/home/runner" ]; then
    chown -R nonroot:nonroot /home/runner 2>/dev/null || true
  fi
  
  if [ -d "/home/agent" ]; then
    chown -R nonroot:nonroot /home/agent 2>/dev/null || true
    chmod -R u+rw /home/agent 2>/dev/null || true
  fi
  
  if [ -d "/home/logs" ]; then
    chown nonroot:nonroot /home/logs 2>/dev/null || true
  fi
  
  if [ -d "/home/submission" ]; then
    chown nonroot:nonroot /home/submission 2>/dev/null || true
  fi
  
  if [ -d "/home/code" ]; then
    chown -R nonroot:nonroot /home/code 2>/dev/null || true
    chmod -R u+rw /home/code 2>/dev/null || true
  fi
  
  # Make sure nonroot can work in /home
  chmod g+w /home 2>/dev/null || true
  
  echo "Minimal permission setup completed at $(date)"
  
  ls -l /home

  echo "Starting grading server at $(date)"
  # Launch grading server, stays alive throughout container lifetime to service agent requests.
  /opt/conda/bin/conda run -n biomlb python /private/grading_server.py
} 2>&1 | tee $LOGS_DIR/entrypoint.log
EOF
        
        chmod +x /home/runner/fast-entrypoint.sh
        
        # Create fast container config with volume mount for fast entrypoint
        mkdir -p /home/runner/biomlbench/environment/config/container_configs/
        cat > /home/runner/biomlbench/environment/config/container_configs/fast.json << 'EOF'
{
    "mem_limit": null,
    "shm_size": "4G",
    "nano_cpus": 16000000000,
    "gpus": 0,
    "fast_entrypoint": "/home/runner/fast-entrypoint.sh"
}
EOF
        """
        
        ssh_setup_cmd = [
            "ssh", "-o", "StrictHostKeyChecking=no",
            f"runner@{public_ip}",
            setup_cmd
        ]
        subprocess.run(ssh_setup_cmd, check=True)
        
        # Optimized prewarming with single finds and better Docker warming
        prewarm_cmd = """
        echo "Starting optimized pre-warming at $(date)"
        
        # Install GNU parallel if not available (much faster than xargs)
        if ! command -v parallel >/dev/null 2>&1; then
            echo "Installing GNU parallel for faster warming..."
            sudo apt-get update -qq && sudo apt-get install -y -qq parallel >/dev/null 2>&1 || echo "Failed to install parallel, using xargs"
        fi
        
        # Use GNU parallel if available for much faster warming
        if command -v parallel >/dev/null 2>&1; then
            echo "Using GNU parallel for faster warming"
            PARALLEL_CMD="parallel -j16 --pipe -N1000 'xargs -0 cat > /dev/null 2>&1'"
        else
            echo "GNU parallel not found, using xargs (still parallel with -P16)"
            PARALLEL_CMD="xargs -0 -P16 -n500 cat > /dev/null 2>&1"
        fi
        
        # Function for fast single-pass warming with progress
        warm_files_fast() {
            local name=$1
            local path=$2
            shift 2
            local find_args="$@"
            
            echo "Pre-warming $name..."
            
            # Single find that streams to parallel processes and shows progress
            find "$path" $find_args -type f -print0 2>/dev/null | \
                tee >(eval $PARALLEL_CMD) | \
                tr '\\0' '\\n' | \
                awk 'BEGIN {count=0} 
                     {count++; if (count % 500 == 0) print "  Progress:", count, "files..."} 
                     END {print "  Found and warmed", count, "files"}'
        }
        
        # 1. Pre-warm Python source files
        warm_files_fast "Python source" "/home/runner/biomlbench/biomlbench/" -name "*.py"
        
        # 2. Pre-warm data cache directories (most critical)
        for dataset in polarishub manual proteingym-dms; do
            if [ -d "/home/runner/.cache/bioml-bench/data/$dataset" ]; then
                warm_files_fast "data cache ($dataset)" "/home/runner/.cache/bioml-bench/data/$dataset"
            fi
        done
        
        # 3. Pre-warm Docker images more efficiently
        echo "Pre-warming Docker images..."
        if [ -d "/var/lib/docker" ]; then
            # Warm only critical Docker metadata and small layers
            {
                # JSON metadata files (instant startup info)
                sudo find /var/lib/docker -name "*.json" -size -1M -print0 2>/dev/null | \
                    xargs -0 -P16 -n100 cat > /dev/null 2>&1
                
                # Small layer files that affect startup
                sudo find /var/lib/docker/overlay2 -path "*/diff/*" -type f -size -5M -print0 2>/dev/null | \
                    head -z -n 1000 | \
                    xargs -0 -P16 -n50 head -c 1M > /dev/null 2>&1
            } &
            DOCKER_PID=$!
        fi
        
        # 4. Pre-warm Python bytecode cache
        warm_files_fast "bytecode cache" "/home/runner/biomlbench" -name "*.pyc"
        
        # 5. Pre-warm key system libraries used by biomlbench
        echo "Pre-warming system libraries..."
        {
            # Python libraries
            find /home/runner/.venv/lib -name "*.so" -size -10M -print0 2>/dev/null | \
                xargs -0 -P16 -n100 head -c 1M > /dev/null 2>&1
        } &
        
        # Wait for background jobs
        wait
        
        echo "Optimized pre-warming completed at $(date)"
        """
        
        # Run SSH with real-time output streaming
        ssh_cmd = [
            "ssh",
            "-o", "StrictHostKeyChecking=no",
            f"runner@{public_ip}",
            prewarm_cmd
        ]
        
        try:
            # Use subprocess.Popen for real-time output
            process = subprocess.Popen(
                ssh_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1  # Line buffered
            )
            
            # Stream output in real-time
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                # Prefix each line with instance ID for clarity
                print(f"[{instance_id}] {line.strip()}")
            
            exit_code = process.wait()
            
            if exit_code == 0:
                log(f"✅ Pre-warmed {instance_id}")
                return True
            else:
                log(f"❌ Failed to pre-warm {instance_id}")
                return False
                
        except Exception as e:
            log(f"❌ Error pre-warming {instance_id}: {str(e)}")
            return False
    
    with ThreadPoolExecutor(max_workers=len(instance_ids)) as executor:
        futures = [executor.submit(warm_instance, instance_id) for instance_id in instance_ids]
        
        # Wait for all futures to complete
        completed = 0
        for future in as_completed(futures):
            try:
                result = future.result()
                completed += 1
                if completed % 5 == 0:  # Progress every 5 instances
                    log(f"Prewarming progress: {completed}/{len(instance_ids)} instances completed")
            except Exception as e:
                log(f"❌ Prewarming failed for one instance: {e}")
                completed += 1
        
        log(f"All {len(instance_ids)} instances processed")

def run_pool_mode(jobs: List[Tuple[str, str]], instance_ids: List[str]):
    """Run jobs using a pool of persistent instances."""
    job_queue = queue.Queue()
    for job in jobs:
        job_queue.put(job)
    
    successful_jobs = 0
    failed_jobs = 0
    lock = threading.Lock()
    
    def worker(instance_id: str):
        nonlocal successful_jobs, failed_jobs
        
        while True:
            try:
                job = job_queue.get(timeout=5)  # Longer timeout for better stability
                agent, task_id = job
                
                log(f"Instance {instance_id} processing: {agent} -> {task_id}")
                success = run_biomlbench_job(instance_id, agent, task_id)
                
                with lock:
                    if success:
                        successful_jobs += 1
                        log(f"🎉 SUCCESS: {agent} -> {task_id}")
                    else:
                        failed_jobs += 1
                        log(f"💥 FAILED: {agent} -> {task_id}")
                
                job_queue.task_done()
                
            except queue.Empty:
                # No more jobs available
                log(f"Instance {instance_id} finished - no more jobs")
                break
    
    # Start worker threads
    threads = []
    for instance_id in instance_ids:
        t = threading.Thread(target=worker, args=(instance_id,))
        t.start()
        threads.append(t)
    
    # Wait for all jobs to complete
    for t in threads:
        t.join()
    
    return successful_jobs, failed_jobs

def main():
    parser = argparse.ArgumentParser(
        description="Deploy BioML-bench jobs on AWS EC2 instances",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Create new instance pool and run jobs (get AMI/SG from setup-aws-resources.sh):
  python scripts/aws-deploy/deploy-aws.py --jobs aws-jobs.txt --s3-prefix s3://biomlbench/v3/artifacts --ami ami-xxxxxxxxx --security-group sg-xxxxxxxxx --concurrent 16
  
  # Use existing instances:
  python scripts/aws-deploy/deploy-aws.py --jobs aws-jobs.txt --s3-prefix s3://my-bucket/experiments/v1 --ami ami-xxxxxxxxx --security-group sg-xxxxxxxxx --existing-instances i-123 i-456 i-789
  
  # GPU instances:
  python scripts/aws-deploy/deploy-aws.py --jobs aws-jobs.txt --s3-prefix s3://my-bucket/data --ami ami-xxxxxxxxx --security-group sg-xxxxxxxxx --concurrent 8 --machine-type gpu
  
Jobs file format (one per line):
  agent,task_id
  aide,polarishub/tdcommons-caco2-wang
  biomni,proteingym-dms/A0A1I9GEU1_NEIME

AWS Resources:
  --ami and --security-group are REQUIRED (no defaults)
  Run setup-aws-resources.sh first to create required resources
  Resource IDs are saved to aws-deploy-config.txt for convenience

This script uses persistent instances that are pre-warmed for optimal performance.
Instances are NOT automatically terminated - clean them up manually when done.
        """
    )
    
    parser.add_argument(
        "--jobs", 
        required=True, 
        help="Path to jobs file (agent,task_id per line)"
    )
    parser.add_argument(
        "--concurrent", 
        type=int, 
        default=15, 
        help="Maximum concurrent instances (default: 15)"
    )
    parser.add_argument(
        "--machine-type", 
        choices=["cpu", "gpu"],
        default="cpu", 
        help="Machine type: cpu (m5.4xlarge, default) or gpu (g4dn.4xlarge)"
    )
    parser.add_argument(
        "--ami",
        required=True,
        help="AMI ID for biomlbench image (create with setup-aws-resources.sh first)"
    )
    parser.add_argument(
        "--key-name",
        default=DEFAULT_KEY_NAME,
        help=f"EC2 key pair name (default: {DEFAULT_KEY_NAME})"
    )
    parser.add_argument(
        "--security-group",
        required=True,
        help="Security group ID (get from setup-aws-resources.sh output or aws-deploy-config.txt)"
    )

    parser.add_argument(
        "--region",
        default=DEFAULT_REGION,
        help=f"AWS region (default: {DEFAULT_REGION})"
    )
    parser.add_argument(
        "--dry-run", 
        action="store_true", 
        help="Show jobs that would be run without executing"
    )
    parser.add_argument(
        "--existing-instances",
        nargs="+",
        help="Instance IDs to use (if not provided, will create new pool)"
    )
    parser.add_argument(
        "--skip-prewarm",
        action="store_true",
        help="Skip pre-warming (useful for already-warmed instances)"
    )
    parser.add_argument(
        "--randomize",
        action="store_true",
        help="Randomize job order before processing"
    )
    parser.add_argument(
        "--tag",
        type=str,
        help="Custom tag to add to instances (e.g., 'run2', 'experiment1')"
    )
    
    args = parser.parse_args()
    
    # Set AWS region
    boto3.setup_default_session(region_name=args.region)
    global ec2
    ec2 = boto3.client('ec2', region_name=args.region)
    
    # Verify AWS setup
    if not args.dry_run:
        if not verify_aws_setup(args.ami, args.key_name, args.security_group):
            log("❌ AWS setup verification failed. Please check your resources.")
            return 1
    
    # Load jobs
    jobs = load_jobs(args.jobs)
    if not jobs:
        log("❌ No valid jobs found")
        return 1
    
    # Randomize jobs if requested
    if args.randomize:
        log("🎲 Randomizing job order...")
        random.shuffle(jobs)
    
    # Dry run
    if args.dry_run:
        log("Dry run - would execute these jobs:")
        for i, (agent, task_id) in enumerate(jobs, 1):
            print(f"  {i}. {agent} -> {task_id}")
        log(f"Instance type: {get_instance_type(args.machine_type)}")
        if args.randomize:
            log("(Job order has been randomized)")
        return 0
    
    # Run jobs using persistent instance pool
    log(f"Machine type: {args.machine_type} ({get_instance_type(args.machine_type)})")
    log(f"Jobs to process: {len(jobs)}")
    
    # Get or create instance pool
    if args.existing_instances:
        instance_ids = args.existing_instances
        log(f"Using existing instances: {instance_ids}")
    else:
        instance_ids = create_instance_pool(
            args.concurrent, args.machine_type, args.ami,
            args.key_name, args.security_group, args.tag
        )
    
    # Pre-warm instances (unless skipped)
    if args.skip_prewarm:
        log("Skipping pre-warming (--skip-prewarm specified)")
    else:
        prewarm_instances(instance_ids)
    
    # Run jobs using pool
    successful_jobs, failed_jobs = run_pool_mode(jobs, instance_ids)
    
    # Summary
    total_jobs = successful_jobs + failed_jobs
    log("=" * 50)
    log("DEPLOYMENT COMPLETE")
    log(f"Total jobs: {total_jobs}")
    log(f"Successful: {successful_jobs}")
    log(f"Failed: {failed_jobs}")
    log(f"Success rate: {successful_jobs/total_jobs*100:.1f}%" if total_jobs > 0 else "N/A")
    
    if successful_jobs > 0:
        log("Check S3 for results: aws s3 ls s3://biomlbench/v3/artifacts/runs/ --recursive")
        log("Check S3 for grades: aws s3 ls s3://biomlbench/v3/artifacts/grades/ --recursive")
    
    return 0 if failed_jobs == 0 else 1

if __name__ == "__main__":
    exit(main()) 