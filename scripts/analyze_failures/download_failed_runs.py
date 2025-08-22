#!/usr/bin/env python3
"""
Download Failed Runs from S3

Script to download and extract all failed run artifacts from biomlbench S3 bucket.
Downloads from s3://biomlbench/v1/artifacts/failed_runs/ and extracts the tar.gz files
to a local directory for analysis.
"""

import argparse
import os
import subprocess
import time
import shutil
from pathlib import Path
from typing import List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

def log(message: str):
    """Simple logging with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def run_command(cmd: List[str], description: str = "", timeout: int = None, cwd: str = None) -> Tuple[int, str]:
    """Run a command and return (exit_code, output)."""
    try:
        if description:
            log(f"Running: {description}")
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            timeout=timeout,
            cwd=cwd
        )
        return result.returncode, result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return 124, "Command timed out"
    except Exception as e:
        return 1, str(e)

def get_failed_run_files(bucket_path: str = "s3://biomlbench/v1/artifacts/failed_runs/") -> List[str]:
    """Get list of all failed run files in S3."""
    log(f"Listing files in {bucket_path}")
    
    exit_code, output = run_command([
        "aws", "s3", "ls", bucket_path, "--recursive"
    ], f"List files in {bucket_path}")
    
    if exit_code != 0:
        log(f"❌ Failed to list S3 files: {output}")
        return []
    
    files = []
    for line in output.strip().split('\n'):
        if not line.strip():
            continue
        # Parse AWS S3 ls output: date time size filepath
        parts = line.split()
        if len(parts) >= 4:
            # Join all parts from index 3 onwards to handle spaces in filenames
            filepath = ' '.join(parts[3:])
            if filepath.endswith('.tar.gz'):
                files.append(f"s3://biomlbench/{filepath}")
    
    log(f"Found {len(files)} failed run files")
    return files

def download_and_extract(s3_path: str, local_dir: Path, keep_compressed: bool = False) -> bool:
    """Download a file from S3 and extract it if it's a tar.gz."""
    filename = Path(s3_path).name
    local_file = local_dir / filename
    
    log(f"Downloading {filename}")
    
    # Download file
    exit_code, output = run_command([
        "aws", "s3", "cp", s3_path, str(local_file)
    ], f"Download {filename}")
    
    if exit_code != 0:
        log(f"❌ Failed to download {filename}: {output}")
        return False
    
    # Extract if it's a tar.gz file
    if filename.endswith('.tar.gz'):
        extract_dir = local_dir / filename.replace('.tar.gz', '')
        extract_dir.mkdir(exist_ok=True)
        
        log(f"Extracting {filename}")
        exit_code, output = run_command([
            "tar", "-xzf", str(local_file), "-C", str(extract_dir)
        ], f"Extract {filename}")
        
        if exit_code != 0:
            log(f"❌ Failed to extract {filename}: {output}")
            return False
        
        log(f"✅ Extracted {filename} to {extract_dir.name}")
        
        # Remove compressed file unless keeping it
        if not keep_compressed:
            local_file.unlink()
            log(f"Removed compressed file {filename}")
    
    return True

def process_failed_run(args: Tuple[str, Path, bool]) -> bool:
    """Process a single failed run (for parallel execution)."""
    s3_path, local_dir, keep_compressed = args
    return download_and_extract(s3_path, local_dir, keep_compressed)

def main():
    parser = argparse.ArgumentParser(
        description="Download and extract failed run artifacts from S3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Download to default directory (./failed_runs):
  python scripts/analyze_failures/download_failed_runs.py
  
  # Download to custom directory:
  python scripts/analyze_failures/download_failed_runs.py --output-dir /path/to/analysis
  
  # Keep compressed files after extraction:
  python scripts/analyze_failures/download_failed_runs.py --keep-compressed
  
  # Use more parallel downloads:
  python scripts/analyze_failures/download_failed_runs.py --concurrent 10
        """
    )
    
    parser.add_argument(
        "--output-dir",
        default="./failed_runs",
        help="Directory to download and extract files to (default: ./failed_runs)"
    )
    parser.add_argument(
        "--bucket-path",
        default="s3://biomlbench/v1/artifacts/failed_runs/",
        help="S3 bucket path to download from (default: s3://biomlbench/v1/artifacts/failed_runs/)"
    )
    parser.add_argument(
        "--concurrent",
        type=int,
        default=5,
        help="Number of concurrent downloads (default: 5)"
    )
    parser.add_argument(
        "--keep-compressed",
        action="store_true",
        help="Keep compressed files after extraction"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without actually doing it"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
        log(f"Output directory: {output_dir.absolute()}")
    
    # Get list of failed run files
    failed_run_files = get_failed_run_files(args.bucket_path)
    if not failed_run_files:
        log("❌ No failed run files found")
        return 1
    
    # Dry run
    if args.dry_run:
        log("Dry run - would download and extract these files:")
        for i, s3_path in enumerate(failed_run_files, 1):
            filename = Path(s3_path).name
            print(f"  {i}. {filename}")
        log(f"Total: {len(failed_run_files)} files")
        return 0
    
    # Download and extract files
    log(f"Starting download with {args.concurrent} concurrent processes")
    log(f"Files to process: {len(failed_run_files)}")
    
    successful_downloads = 0
    failed_downloads = 0
    
    # Prepare arguments for parallel processing
    process_args = [
        (s3_path, output_dir, args.keep_compressed) 
        for s3_path in failed_run_files
    ]
    
    with ThreadPoolExecutor(max_workers=args.concurrent) as executor:
        # Submit all download tasks
        future_to_path = {
            executor.submit(process_failed_run, arg): arg[0] 
            for arg in process_args
        }
        
        # Process completed downloads
        for future in as_completed(future_to_path):
            s3_path = future_to_path[future]
            filename = Path(s3_path).name
            
            try:
                success = future.result()
                if success:
                    successful_downloads += 1
                    log(f"🎉 SUCCESS: {filename}")
                else:
                    failed_downloads += 1
                    log(f"💥 FAILED: {filename}")
            except Exception as e:
                failed_downloads += 1
                log(f"💥 EXCEPTION: {filename}: {e}")
    
    # Summary
    total_files = successful_downloads + failed_downloads
    log("=" * 50)
    log("DOWNLOAD COMPLETE")
    log(f"Total files: {total_files}")
    log(f"Successful: {successful_downloads}")
    log(f"Failed: {failed_downloads}")
    log(f"Success rate: {successful_downloads/total_files*100:.1f}%" if total_files > 0 else "N/A")
    
    if successful_downloads > 0:
        log(f"Failed runs extracted to: {output_dir.absolute()}")
        log("You can now analyze the failed runs in the extracted directories")
        
        # Show directory structure
        if output_dir.exists():
            log("\nDirectory structure:")
            for item in sorted(output_dir.iterdir()):
                if item.is_dir():
                    file_count = len(list(item.rglob('*')))
                    log(f"  📁 {item.name}/ ({file_count} files)")
    
    return 0 if failed_downloads == 0 else 1

if __name__ == "__main__":
    exit(main()) 