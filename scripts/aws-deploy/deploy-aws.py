#!/usr/bin/env python3
"""
AWS Deployment Script for BioML-bench

AWS equivalent of the GCP deployment system that:
1. Reads jobs from a text file (can be run from project root)
2. Creates EC2 instances with biomlbench AMI
3. Runs biomlbench pipeline on each instance
4. Verifies S3 uploads
5. Cleans up instances
6. Manages parallel execution
"""

import argparse
import subprocess
import time
import uuid
import json
import boto3
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple, Optional

# AWS clients
ec2 = boto3.client('ec2')
ssm = boto3.client('ssm')

# Default configuration
DEFAULT_REGION = "us-east-1"
DEFAULT_AMI = 'ami-0cbfb39216f82ad5f'  # Will be set from environment or argument
DEFAULT_KEY_NAME = "biomlbench-key"
DEFAULT_SECURITY_GROUP = "sg-050d87e9eef71d1ce"

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
                              max_retries: int = None) -> Optional[str]:
    """Create an EC2 instance, retrying until success. Returns instance ID if successful."""
    log(f"Creating instance: {instance_name} (machine type: {machine_type})")
    
    instance_type = get_instance_type(machine_type)
    
    # User data script to set up the instance
    user_data = """#!/bin/bash
# Set up SSH for the runner user
mkdir -p /home/runner/.ssh
echo 'ssh-rsa YOUR_PUBLIC_KEY_HERE' >> /home/runner/.ssh/authorized_keys
chown -R runner:runner /home/runner/.ssh
chmod 700 /home/runner/.ssh
chmod 600 /home/runner/.ssh/authorized_keys

# Ensure AWS CLI is configured (if not already in AMI)
aws configure set region us-east-1
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
                        ]
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

def wait_for_ssm(instance_id: str, timeout: int = 300) -> bool:
    """Wait for SSM agent to become available on the instance."""
    log(f"Waiting for SSM on instance {instance_id}...")
    
    start_time = time.time()
    while time.time() - start_time < timeout:
        try:
            response = ssm.describe_instance_information(
                Filters=[
                    {'Key': 'InstanceIds', 'Values': [instance_id]}
                ]
            )
            
            if response['InstanceInformationList']:
                if response['InstanceInformationList'][0]['PingStatus'] == 'Online':
                    log(f"✅ SSM ready on instance {instance_id}")
                    return True
        except Exception:
            pass
        
        time.sleep(10)
    
    log(f"❌ SSM timeout on instance {instance_id}")
    return False

def run_biomlbench_job(instance_id: str, agent: str, task_id: str) -> bool:
    """Run the biomlbench pipeline on an instance using SSM."""
    log(f"Running job on {instance_id}: {agent} -> {task_id}")
    
    # Build the command to run on the instance
    remote_cmd = f"""#!/bin/bash
    set -e
    cd /home/runner/biomlbench
    source .venv/bin/activate
    
    # Run agent (CPU-only mode)
    biomlbench run-agent --agent {agent} --task-id {task_id} --cpu-only
    
    # Get the specific run group ID that was just created
    LATEST_RUN=$(find runs/ -name "*run-group_{agent}" -type d | sort | tail -1)
    RUN_GROUP_ID=$(basename "$LATEST_RUN")
    echo "📁 Run group: $RUN_GROUP_ID"
    
    # Grade results
    biomlbench grade --submission "$LATEST_RUN/submission.jsonl" --output-dir results/
    
    # Get the grading timestamp from the most recent grading report
    GRADING_REPORT=$(find results/ -name "*_grading_report.json" | sort | tail -1)
    GRADING_TIMESTAMP=$(basename "$GRADING_REPORT" | cut -d'_' -f1)
    echo "📊 Grading timestamp: $GRADING_TIMESTAMP"
    
    # Convert task_id to S3-safe format for the organized structure
    TASK_ID_SAFE=$(echo "{task_id}" | sed 's/\//-/g' | sed 's/_/-/g')
    
    # Show the exact S3 paths for this specific run
    echo "📤 S3 artifacts for this run:"
    echo "  Run artifacts:"
    echo "    s3://biomlbench/v2/artifacts/runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz"
    echo "    OR s3://biomlbench/v2/artifacts/failed_runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz (if failed)"
    echo "  Grading artifacts:"
    echo "    s3://biomlbench/v2/artifacts/grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_grading_report.json.gz"
    echo "    s3://biomlbench/v2/artifacts/grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_individual_reports.tar.gz"
    echo "    OR s3://biomlbench/v2/artifacts/failed_grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_*.gz (if failed)"
    
    # Verify these specific paths exist
    echo "🔍 Verifying uploads..."
    if aws s3 ls s3://biomlbench/v2/artifacts/runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Run artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v2/artifacts/runs/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Run artifacts uploaded successfully (flat structure)"
    elif aws s3 ls s3://biomlbench/v2/artifacts/failed_runs/{agent}/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Failed run artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v2/artifacts/failed_runs/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Failed run artifacts uploaded successfully (flat structure)"
    else
        echo "❌ No run artifacts found in S3!"
        exit 1
    fi
    
    if aws s3 ls s3://biomlbench/v2/artifacts/grades/{agent}/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_grading_report.json.gz > /dev/null 2>&1; then
        echo "✅ Grading artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v2/artifacts/grades/${{GRADING_TIMESTAMP}}_grading_report.json.gz > /dev/null 2>&1; then
        echo "✅ Grading artifacts uploaded successfully (flat structure)"
    elif aws s3 ls s3://biomlbench/v2/artifacts/failed_grades/{agent}/$TASK_ID_SAFE/ | grep -q "$GRADING_TIMESTAMP" > /dev/null 2>&1; then
        echo "✅ Failed grading artifacts uploaded successfully (organized structure)"
    elif aws s3 ls s3://biomlbench/v2/artifacts/failed_grades/ | grep -q "$GRADING_TIMESTAMP" > /dev/null 2>&1; then
        echo "✅ Failed grading artifacts uploaded successfully (flat structure)"
    else
        echo "❌ No grading artifacts found in S3!"
        exit 1
    fi
    
    echo "✅ Job completed successfully"
    """
    
    try:
        # Send command via SSM
        response = ssm.send_command(
            InstanceIds=[instance_id],
            DocumentName="AWS-RunShellScript",
            Parameters={
                'commands': [remote_cmd],
                'workingDirectory': ['/home/runner'],
                'executionTimeout': ['3600']  # 1 hour timeout
            }
        )
        
        command_id = response['Command']['CommandId']
        
        # Wait for command to complete
        time.sleep(5)  # Give it a moment to start
        
        while True:
            try:
                output = ssm.get_command_invocation(
                    CommandId=command_id,
                    InstanceId=instance_id
                )
                
                if output['Status'] in ['Success', 'Failed', 'Cancelled', 'TimedOut']:
                    if output['Status'] == 'Success':
                        log(f"✅ Job completed on {instance_id}: {agent} -> {task_id}")
                        # Show the output
                        if output['StandardOutputContent']:
                            print("\n" + "="*60)
                            print("JOB OUTPUT:")
                            print("="*60)
                            print(output['StandardOutputContent'])
                            print("="*60 + "\n")
                        return True
                    else:
                        log(f"❌ Job failed on {instance_id}: {agent} -> {task_id}")
                        log(f"Status: {output['Status']}")
                        if output['StandardErrorContent']:
                            log(f"Error: {output['StandardErrorContent']}")
                        return False
                
                time.sleep(10)
                
            except ssm.exceptions.InvocationDoesNotExist:
                time.sleep(10)
                
    except Exception as e:
        log(f"❌ Failed to run job on {instance_id}: {str(e)}")
        return False

def delete_instance(instance_id: str):
    """Delete an EC2 instance."""
    log(f"Deleting instance: {instance_id}")
    
    try:
        ec2.terminate_instances(InstanceIds=[instance_id])
        log(f"✅ Deleted instance: {instance_id}")
    except Exception as e:
        log(f"⚠️  Failed to delete instance {instance_id}: {str(e)}")

def process_job(job: Tuple[str, str], machine_type: str, ami_id: str,
                key_name: str = DEFAULT_KEY_NAME,
                security_group: str = DEFAULT_SECURITY_GROUP) -> bool:
    """Process a single job (agent, task_id) on a dedicated instance."""
    agent, task_id = job
    
    # Generate unique instance name (EC2 max 255 chars for tags)
    safe_task = task_id.replace('/', '-').replace('_', '-').lower()
    uuid_suffix = uuid.uuid4().hex[:8]
    
    # Calculate max length for task part
    max_task_len = 255 - 7 - len(agent) - 2 - 8  # "bioml-" + agent + "-" + "-" + uuid
    if len(safe_task) > max_task_len:
        # Truncate and add hash for uniqueness
        import hashlib
        task_hash = hashlib.md5(safe_task.encode()).hexdigest()[:4]
        safe_task = safe_task[:max_task_len-5] + '-' + task_hash
    
    instance_name = f"bioml-{agent}-{safe_task}-{uuid_suffix}"
    instance_id = None
    
    try:
        # Create instance (max 3 attempts)
        instance_id = create_instance_with_retry(
            instance_name, machine_type, ami_id,
            key_name, security_group,
            max_retries=3
        )
        if not instance_id:
            return False
        
        # Wait for SSM
        if not wait_for_ssm(instance_id):
            return False
        
        # Run job
        success = run_biomlbench_job(instance_id, agent, task_id)
        
        return success
        
    finally:
        # Always clean up instance
        if instance_id:
            delete_instance(instance_id)

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

def main():
    parser = argparse.ArgumentParser(
        description="Deploy BioML-bench jobs on AWS EC2 instances",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # From project root with CPU (default):
  python scripts/aws-deploy/deploy-aws.py --jobs aws-jobs.txt --concurrent 5 --ami ami-xxxxxxxxx
  
  # With GPU machines:
  python scripts/aws-deploy/deploy-aws.py --jobs aws-jobs.txt --concurrent 5 --machine-type gpu --ami ami-xxxxxxxxx
  
  # From scripts/aws-deploy directory:
  python deploy-aws.py --jobs jobs.txt --concurrent 5 --ami ami-xxxxxxxxx
  
Jobs file format (one per line):
  agent,task_id
  aide,polarishub/tdcommons-caco2-wang
  biomni,proteingym-dms/A0A1I9GEU1_NEIME
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
        default=DEFAULT_AMI,
        help="AMI ID for biomlbench image (required)"
    )
    parser.add_argument(
        "--key-name",
        default=DEFAULT_KEY_NAME,
        help=f"EC2 key pair name (default: {DEFAULT_KEY_NAME})"
    )
    parser.add_argument(
        "--security-group",
        default=DEFAULT_SECURITY_GROUP,
        help=f"Security group ID (default: {DEFAULT_SECURITY_GROUP})"
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
    
    args = parser.parse_args()
    
    # Set AWS region
    boto3.setup_default_session(region_name=args.region)
    global ec2, ssm
    ec2 = boto3.client('ec2', region_name=args.region)
    ssm = boto3.client('ssm', region_name=args.region)
    
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
    
    # Dry run
    if args.dry_run:
        log("Dry run - would execute these jobs:")
        for i, (agent, task_id) in enumerate(jobs, 1):
            print(f"  {i}. {agent} -> {task_id}")
        log(f"Instance type: {get_instance_type(args.machine_type)}")
        return 0
    
    # Run jobs
    log(f"Starting deployment with {args.concurrent} concurrent instances")
    log(f"Machine type: {args.machine_type} ({get_instance_type(args.machine_type)})")
    log(f"Jobs to process: {len(jobs)}")
    
    successful_jobs = 0
    failed_jobs = 0
    
    with ThreadPoolExecutor(max_workers=args.concurrent) as executor:
        # Submit all jobs
        future_to_job = {
            executor.submit(
                process_job, job, args.machine_type, args.ami,
                args.key_name, args.security_group
            ): job 
            for job in jobs
        }
        
        # Process completed jobs
        for future in as_completed(future_to_job):
            job = future_to_job[future]
            agent, task_id = job
            
            try:
                success = future.result()
                if success:
                    successful_jobs += 1
                    log(f"🎉 SUCCESS: {agent} -> {task_id}")
                else:
                    failed_jobs += 1
                    log(f"💥 FAILED: {agent} -> {task_id}")
            except Exception as e:
                failed_jobs += 1
                log(f"💥 EXCEPTION: {agent} -> {task_id}: {e}")
    
    # Summary
    total_jobs = successful_jobs + failed_jobs
    log("=" * 50)
    log("DEPLOYMENT COMPLETE")
    log(f"Total jobs: {total_jobs}")
    log(f"Successful: {successful_jobs}")
    log(f"Failed: {failed_jobs}")
    log(f"Success rate: {successful_jobs/total_jobs*100:.1f}%" if total_jobs > 0 else "N/A")
    
    if successful_jobs > 0:
        log("Check S3 for results: aws s3 ls s3://biomlbench/v2/artifacts/runs/ --recursive")
        log("Check S3 for grades: aws s3 ls s3://biomlbench/v2/artifacts/grades/ --recursive")
    
    return 0 if failed_jobs == 0 else 1

if __name__ == "__main__":
    exit(main()) 