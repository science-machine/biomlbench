#!/usr/bin/env python3
"""
GCP Deployment Script for BioML-bench

Simple orchestrator that:
1. Reads jobs from a text file (can be run from project root)
2. Creates GCP VMs with biomlbench image
3. Runs biomlbench pipeline on each VM
4. Verifies S3 uploads
5. Cleans up VMs
6. Manages parallel execution
"""

import argparse
import subprocess
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

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

def create_vm_with_retry(vm_name: str, zone: str = "us-central1-a", max_retries: int = None) -> bool:
    """Create a VM, retrying until success."""
    log(f"Creating VM: {vm_name}")
    
    cmd = [
        "gcloud", "compute", "instances", "create", vm_name,
        "--zone", zone,
        "--machine-type", "g2-standard-12",
        "--maintenance-policy", "TERMINATE", 
        "--image", "biomlbenchhm",
        "--boot-disk-size", "500G"
    ]
    
    attempt = 0
    while True:
        attempt += 1
        if max_retries and attempt > max_retries:
            log(f"❌ Failed to create {vm_name} after {max_retries} attempts")
            return False
            
        exit_code, output = run_command(cmd, f"Create VM {vm_name} (attempt {attempt})")
        
        if exit_code == 0:
            log(f"✅ Created VM: {vm_name}")
            return True
        else:
            log(f"⚠️  VM creation failed (attempt {attempt}): {output.strip()}")
            time.sleep(5)  # Wait before retry

def wait_for_ssh(vm_name: str, zone: str = "us-central1-a", timeout: int = 300) -> bool:
    """Wait for SSH to become available on the VM."""
    log(f"Waiting for SSH on {vm_name}...")
    
    start_time = time.time()
    while time.time() - start_time < timeout:
        exit_code, _ = run_command([
            "gcloud", "compute", "ssh", f"runner@{vm_name}",
            "--zone", zone,
            "--command", "echo 'SSH ready'",
            "--quiet"
        ], timeout=10)
        
        if exit_code == 0:
            log(f"✅ SSH ready on {vm_name}")
            return True
        
        time.sleep(10)
    
    log(f"❌ SSH timeout on {vm_name}")
    return False

def run_biomlbench_job(vm_name: str, agent: str, task_id: str, zone: str = "us-central1-a") -> bool:
    """Run the biomlbench pipeline on a VM."""
    log(f"Running job on {vm_name}: {agent} -> {task_id}")
    
    # Build the command to run on the VM
    remote_cmd = f"""
    set -e
    cd biomlbench
    source .venv/bin/activate
    
    # Run agent
    biomlbench run-agent --agent {agent} --task-id {task_id}
    
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
    TASK_ID_SAFE=$(echo "$task_id" | sed 's/\//-/g' | sed 's/_/-/g')
    
    # Show the exact S3 paths for this specific run
    echo "📤 S3 artifacts for this run:"
    echo "  Run artifacts:"
    echo "    s3://biomlbench/v1/artifacts/runs/$agent/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz"
    echo "    OR s3://biomlbench/v1/artifacts/failed_runs/$agent/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz (if failed)"
    echo "  Grading artifacts:"
    echo "    s3://biomlbench/v1/artifacts/grades/$agent/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_grading_report.json.gz"
    echo "    s3://biomlbench/v1/artifacts/grades/$agent/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_individual_reports.tar.gz"
    echo "    OR s3://biomlbench/v1/artifacts/failed_grades/$agent/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_*.gz (if failed)"
    
    # Verify these specific paths exist
    echo "🔍 Verifying uploads..."
    if aws s3 ls s3://biomlbench/v1/artifacts/runs/$agent/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Run artifacts uploaded successfully"
    elif aws s3 ls s3://biomlbench/v1/artifacts/failed_runs/$agent/$TASK_ID_SAFE/$RUN_GROUP_ID.tar.gz > /dev/null 2>&1; then
        echo "✅ Failed run artifacts uploaded successfully"
    else
        echo "❌ No run artifacts found in S3!"
        exit 1
    fi
    
    if aws s3 ls s3://biomlbench/v1/artifacts/grades/$agent/$TASK_ID_SAFE/${{GRADING_TIMESTAMP}}_grading_report.json.gz > /dev/null 2>&1; then
        echo "✅ Grading artifacts uploaded successfully"
    elif aws s3 ls s3://biomlbench/v1/artifacts/failed_grades/$agent/$TASK_ID_SAFE/ | grep -q "$GRADING_TIMESTAMP" > /dev/null 2>&1; then
        echo "✅ Failed grading artifacts uploaded successfully"
    else
        echo "❌ No grading artifacts found in S3!"
        exit 1
    fi
    
    echo "✅ Job completed successfully"
    """
    
    exit_code, output = run_command([
        "gcloud", "compute", "ssh", f"runner@{vm_name}",
        "--zone", zone,
        "--command", remote_cmd,
        "--quiet"
    ])
    
    if exit_code == 0:
        log(f"✅ Job completed on {vm_name}: {agent} -> {task_id}")
        # Show the output so we can see S3 paths
        if output.strip():
            print("\n" + "="*60)
            print("JOB OUTPUT:")
            print("="*60)
            print(output.strip())
            print("="*60 + "\n")
        return True
    else:
        log(f"❌ Job failed on {vm_name}: {agent} -> {task_id}")
        log(f"Error output: {output}")
        return False

def delete_vm(vm_name: str, zone: str = "us-central1-a"):
    """Delete a VM."""
    log(f"Deleting VM: {vm_name}")
    
    exit_code, output = run_command([
        "gcloud", "compute", "instances", "delete", vm_name,
        "--zone", zone,
        "--quiet"
    ])
    
    if exit_code == 0:
        log(f"✅ Deleted VM: {vm_name}")
    else:
        log(f"⚠️  Failed to delete VM {vm_name}: {output}")

def process_job(job: Tuple[str, str], zone: str = "us-central1-a") -> bool:
    """Process a single job (agent, task_id) on a dedicated VM."""
    agent, task_id = job
    
    # Generate unique VM name (GCP max 63 chars, lowercase, no underscores)
    safe_task = task_id.replace('/', '-').replace('_', '-').lower()
    uuid_suffix = uuid.uuid4().hex[:8]
    
    # Calculate max length for task part: 63 - "bioml-" - agent - "-" - "-" - uuid = 63 - 6 - len(agent) - 2 - 8
    max_task_len = 63 - 6 - len(agent) - 2 - 8
    if len(safe_task) > max_task_len:
        # Truncate and add hash for uniqueness
        import hashlib
        task_hash = hashlib.md5(safe_task.encode()).hexdigest()[:4]
        safe_task = safe_task[:max_task_len-5] + '-' + task_hash
    
    vm_name = f"bioml-{agent}-{safe_task}-{uuid_suffix}"
    
    try:
        # Create VM (max 3 attempts)
        if not create_vm_with_retry(vm_name, zone, max_retries=3):
            return False
        
        # Wait for SSH
        if not wait_for_ssh(vm_name, zone):
            delete_vm(vm_name, zone)
            return False
        
        # Run job
        success = run_biomlbench_job(vm_name, agent, task_id, zone)
        
        return success
        
    finally:
        # Always clean up VM
        delete_vm(vm_name, zone)

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

def main():
    parser = argparse.ArgumentParser(
        description="Deploy BioML-bench jobs on GCP VMs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # From project root:
  python scripts/gcp-deploy/deploy.py --jobs gcp-jobs.txt --concurrent 5
  
  # From scripts/gcp-deploy directory:
  python deploy.py --jobs jobs.txt --concurrent 5
  
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
        help="Maximum concurrent VMs (default: 15, max recommended: 16 based on L4 GPU quota)"
    )
    parser.add_argument(
        "--zone", 
        default="us-central1-a", 
        help="GCP zone (default: us-central1-a)"
    )
    parser.add_argument(
        "--dry-run", 
        action="store_true", 
        help="Show jobs that would be run without executing"
    )
    
    args = parser.parse_args()
    
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
        return 0
    
    # Run jobs
    log(f"Starting deployment with {args.concurrent} concurrent VMs")
    log(f"Jobs to process: {len(jobs)}")
    
    successful_jobs = 0
    failed_jobs = 0
    
    with ThreadPoolExecutor(max_workers=args.concurrent) as executor:
        # Submit all jobs
        future_to_job = {
            executor.submit(process_job, job, args.zone): job 
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
        log("Check S3 for results: aws s3 ls s3://biomlbench/v1/artifacts/runs/ --recursive")
        log("Check S3 for grades: aws s3 ls s3://biomlbench/v1/artifacts/grades/ --recursive")
    
    return 0 if failed_jobs == 0 else 1

if __name__ == "__main__":
    exit(main()) 