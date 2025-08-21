#!/usr/bin/env python3
"""
BioML-bench Dummy Agent Artifact Cleanup Script

Finds and deletes all dummy agent artifacts in s3://biomlbench/v1/artifacts/

Features:
- Dry-run mode to preview what would be deleted
- Interactive confirmation before deletion
- Handles both organized and flat artifact structures
- Parallel deletion for efficiency
- Comprehensive logging

Usage Examples:
  # Dry run to see what would be deleted
  python scripts/cleanup_dummy_artifacts.py --dry-run
  
  # Interactive deletion with confirmation
  python scripts/cleanup_dummy_artifacts.py
  
  # Force deletion without confirmation (BE CAREFUL!)
  python scripts/cleanup_dummy_artifacts.py --force
  
  # Custom bucket/prefix
  python scripts/cleanup_dummy_artifacts.py --bucket my-bucket --prefix v2/artifacts
"""

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Set
import time

def log(message: str):
    """Simple logging with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def run_aws_command(cmd: List[str], description: str = "") -> tuple[int, str]:
    """Run AWS CLI command and return (exit_code, output)."""
    try:
        if description:
            log(f"Running: {description}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.returncode, result.stdout + result.stderr
    except Exception as e:
        return 1, str(e)

def discover_dummy_artifacts(bucket: str, prefix: str) -> Set[str]:
    """Discover all dummy agent artifacts in S3."""
    log("🔍 Discovering dummy agent artifacts...")
    
    dummy_artifacts = set()
    
    # Define artifact types to search
    artifact_types = ['runs', 'grades', 'failed_runs', 'failed_grades']
    
    for artifact_type in artifact_types:
        s3_path = f"s3://{bucket}/{prefix}/{artifact_type}/"
        log(f"  Searching {artifact_type}...")
        
        exit_code, output = run_aws_command([
            "aws", "s3", "ls", s3_path, "--recursive"
        ], f"List {artifact_type}")
        
        if exit_code != 0:
            log(f"⚠️  Failed to list {artifact_type}: {output}")
            continue
        
        for line in output.strip().split('\n'):
            if not line.strip():
                continue
                
            # Parse S3 ls output: "2025-08-20 01:10:23    123 path/to/file"
            parts = line.split()
            if len(parts) < 4:
                continue
                
            s3_key = ' '.join(parts[3:])
            full_path = f"s3://{bucket}/{s3_key}"
            
            # Check if this is a dummy agent artifact
            is_dummy_artifact = False
            
            # Method 1: Organized structure (artifacts/type/agent/task/file)
            key_parts = s3_key.split('/')
            if len(key_parts) >= 4:
                # Check if agent directory is 'dummy'
                if key_parts[-3] == 'dummy':
                    is_dummy_artifact = True
            
            # Method 2: Flat structure (artifacts/type/timestamp_run-group_dummy.ext)
            filename = key_parts[-1]
            if '_run-group_dummy.' in filename or filename.endswith('_dummy.tar.gz') or filename.endswith('_dummy.json.gz'):
                is_dummy_artifact = True
            
            # Method 3: Directory structure contains dummy
            if '/dummy/' in s3_key:
                is_dummy_artifact = True
            
            if is_dummy_artifact:
                dummy_artifacts.add(full_path)
    
    log(f"📊 Found {len(dummy_artifacts)} dummy agent artifacts")
    return dummy_artifacts

def delete_single_artifact(s3_path: str) -> tuple[str, bool, str]:
    """Delete a single artifact from S3."""
    try:
        exit_code, output = run_aws_command([
            "aws", "s3", "rm", s3_path
        ])
        
        success = exit_code == 0
        return s3_path, success, output.strip() if output.strip() else "Success"
        
    except Exception as e:
        return s3_path, False, str(e)

def delete_artifacts_parallel(artifacts: Set[str], max_workers: int = 10) -> tuple[int, int]:
    """Delete artifacts in parallel."""
    log(f"🗑️  Deleting {len(artifacts)} artifacts with {max_workers} parallel workers...")
    
    successful_deletions = 0
    failed_deletions = 0
    completed = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all deletion tasks
        future_to_path = {
            executor.submit(delete_single_artifact, s3_path): s3_path
            for s3_path in artifacts
        }
        
        # Process completed deletions
        for future in as_completed(future_to_path):
            s3_path = future_to_path[future]
            completed += 1
            
            try:
                path, success, message = future.result()
                if success:
                    successful_deletions += 1
                else:
                    failed_deletions += 1
                    log(f"❌ Failed to delete {path}: {message}")
            except Exception as e:
                failed_deletions += 1
                log(f"❌ Exception deleting {s3_path}: {e}")
            
            # Progress report every 25 deletions
            if completed % 25 == 0:
                log(f"Progress: {completed}/{len(artifacts)} deletions completed ({completed/len(artifacts)*100:.1f}%)")
    
    return successful_deletions, failed_deletions

def confirm_deletion(artifacts: Set[str]) -> bool:
    """Ask user to confirm deletion."""
    print("\n" + "="*80)
    print("🚨 DELETION CONFIRMATION")
    print("="*80)
    print(f"You are about to DELETE {len(artifacts)} dummy agent artifacts from S3.")
    print("This operation is PERMANENT and CANNOT be undone!")
    print("\nFirst 10 artifacts that will be deleted:")
    for i, artifact in enumerate(sorted(artifacts)):
        if i >= 10:
            print(f"... and {len(artifacts) - 10} more")
            break
        print(f"  - {artifact}")
    
    print("\n" + "="*80)
    
    while True:
        response = input("Are you sure you want to delete these artifacts? (yes/no): ").strip().lower()
        if response in ['yes', 'y']:
            return True
        elif response in ['no', 'n']:
            return False
        else:
            print("Please enter 'yes' or 'no'")

def main():
    parser = argparse.ArgumentParser(
        description="Delete all dummy agent artifacts from BioML-bench S3 bucket"
    )
    parser.add_argument(
        "--bucket",
        default="biomlbench",
        help="S3 bucket name (default: biomlbench)"
    )
    parser.add_argument(
        "--prefix",
        default="v1/artifacts",
        help="S3 prefix (default: v1/artifacts)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without actually deleting"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Skip confirmation prompt (BE CAREFUL!)"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=10,
        help="Number of parallel deletion workers (default: 10)"
    )
    
    args = parser.parse_args()
    
    # Discover dummy artifacts
    try:
        artifacts = discover_dummy_artifacts(args.bucket, args.prefix)
    except Exception as e:
        log(f"❌ Failed to discover artifacts: {e}")
        return 1
    
    if not artifacts:
        log("✅ No dummy agent artifacts found!")
        return 0
    
    # Group artifacts by type for summary
    artifact_summary = {}
    for artifact in artifacts:
        if '/runs/' in artifact:
            artifact_type = 'runs'
        elif '/grades/' in artifact:
            artifact_type = 'grades'
        elif '/failed_runs/' in artifact:
            artifact_type = 'failed_runs'
        elif '/failed_grades/' in artifact:
            artifact_type = 'failed_grades'
        else:
            artifact_type = 'other'
        
        if artifact_type not in artifact_summary:
            artifact_summary[artifact_type] = 0
        artifact_summary[artifact_type] += 1
    
    # Show summary
    print("\n" + "="*60)
    print("📋 DUMMY AGENT ARTIFACTS SUMMARY")
    print("="*60)
    for artifact_type, count in sorted(artifact_summary.items()):
        print(f"  {artifact_type}: {count} artifacts")
    print(f"  TOTAL: {len(artifacts)} artifacts")
    
    if args.dry_run:
        print("\n🔍 DRY RUN - No artifacts will be deleted")
        print("\nAll dummy agent artifacts that would be deleted:")
        for artifact in sorted(artifacts):
            print(f"  - {artifact}")
        return 0
    
    # Confirm deletion
    if not args.force:
        if not confirm_deletion(artifacts):
            log("❌ Deletion cancelled by user")
            return 1
    
    # Perform deletion
    log("🚀 Starting deletion process...")
    start_time = time.time()
    
    try:
        successful, failed = delete_artifacts_parallel(artifacts, args.max_workers)
        
        elapsed = time.time() - start_time
        log(f"✅ Deletion completed in {elapsed:.1f} seconds")
        log(f"📊 Successfully deleted: {successful}")
        if failed > 0:
            log(f"❌ Failed deletions: {failed}")
            return 1
        else:
            log("🎉 All dummy agent artifacts successfully deleted!")
            return 0
            
    except KeyboardInterrupt:
        log("❌ Deletion interrupted by user")
        return 1
    except Exception as e:
        log(f"❌ Deletion failed: {e}")
        return 1

if __name__ == "__main__":
    exit(main()) 