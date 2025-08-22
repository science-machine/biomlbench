#!/usr/bin/env python3
"""
Parse AIDE Journal Files

Script to extract and format AIDE agent journal files from failed runs for human review.
Finds all journal.json files in the failed_runs directory and converts them to readable
text formats organized by agent, task, and run ID.
"""

import argparse
import json
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

def log(message: str):
    """Simple logging with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def find_aide_journal_files(failed_runs_dir: str = "./failed_runs") -> List[Tuple[str, Dict]]:
    """Find all AIDE journal files in failed runs directory."""
    journal_files = []
    failed_runs_path = Path(failed_runs_dir)
    
    if not failed_runs_path.exists():
        log(f"❌ Failed runs directory not found: {failed_runs_dir}")
        return []
    
    log(f"Searching for AIDE journal files in {failed_runs_path.absolute()}")
    
    # Pattern: failed_runs/<run_id>/<run_id>/<domain>/<task_id>_<uuid>/logs/<unique_name>/journal.json
    for journal_path in failed_runs_path.rglob("**/logs/*/journal.json"):
        # Extract run information from path
        path_parts = journal_path.parts
        
        # Find the run_id (should appear twice in the path)
        run_id = None
        for i, part in enumerate(path_parts):
            if part.startswith("2025-") and "run-group_aide" in part:
                run_id = part
                break
        
        if run_id and "aide" in run_id.lower():
            # Extract additional metadata from path
            try:
                # Get domain and task info from path
                logs_index = path_parts.index("logs")
                if logs_index >= 3:
                    domain = path_parts[logs_index - 2]  # e.g., "kaggle", "manual"
                    task_with_uuid = path_parts[logs_index - 1]  # e.g., "uw-madison-gi-tract-image-segmentation_7e12c4d2..."
                    task_id = task_with_uuid.rsplit('_', 1)[0]  # Remove UUID suffix
                    
                    metadata = {
                        'run_id': run_id,
                        'domain': domain,
                        'task_id': task_id,
                        'task_with_uuid': task_with_uuid,
                        'logs_dir': path_parts[logs_index + 1]
                    }
                    journal_files.append((str(journal_path), metadata))
                    
            except (ValueError, IndexError):
                log(f"⚠️  Could not parse path structure: {journal_path}")
                continue
    
    log(f"Found {len(journal_files)} AIDE journal files")
    return journal_files

def parse_journal_file(journal_path: str) -> Optional[Dict]:
    """Parse a journal.json file."""
    try:
        with open(journal_path, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        log(f"❌ Failed to parse {journal_path}: {e}")
        return None

def format_step_text(node: Dict, step_num: int) -> str:
    """Format a single step as plain text."""
    text_parts = []
    
    # Step header
    ctime_readable = datetime.fromtimestamp(node.get('ctime', 0)).strftime('%Y-%m-%d %H:%M:%S') if node.get('ctime') else 'Unknown'
    text_parts.append("=" * 100)
    text_parts.append(f"STEP {step_num:2d} - {ctime_readable}")
    text_parts.append("=" * 100)
    
    # Plan
    if node.get('plan'):
        text_parts.append(f"\n📋 PLAN:")
        text_parts.append("-" * 50)
        text_parts.append(str(node['plan']))
    
    # Code
    if node.get('code'):
        text_parts.append(f"\n💻 CODE:")
        text_parts.append("-" * 50)
        text_parts.append(str(node['code']))
    
    # Terminal Output
    if node.get('_term_out'):
        term_out = node['_term_out']
        if isinstance(term_out, list):
            term_out = '\n'.join(term_out)
        text_parts.append(f"\n📤 OUTPUT:")
        text_parts.append("-" * 50)
        text_parts.append(str(term_out))
    
    # Analysis
    if node.get('analysis'):
        text_parts.append(f"\n🔍 ANALYSIS:")
        text_parts.append("-" * 50)
        text_parts.append(str(node['analysis']))
    
    # Error information
    if node.get('exc_type') or node.get('exc_info'):
        exc_type = node.get('exc_type', 'Unknown')
        exc_info = node.get('exc_info', 'No details')
        text_parts.append(f"\n❌ ERROR:")
        text_parts.append("-" * 50)
        text_parts.append(f"Type: {exc_type}")
        text_parts.append(f"Info: {exc_info}")
    
    # Execution time and status
    details = []
    if node.get('exec_time') is not None:
        details.append(f"⏱️  Execution time: {node['exec_time']:.2f}s")
    if node.get('is_buggy'):
        details.append("⚠️  Marked as buggy")
    if node.get('metric') is not None:
        details.append(f"📊 Metric: {node['metric']}")
    
    if details:
        text_parts.append(f"\n📊 DETAILS:")
        text_parts.append("-" * 50)
        text_parts.append('\n'.join(details))
    
    return '\n'.join(text_parts)

def process_journal(args: Tuple[str, Dict]) -> bool:
    """Process a single journal file."""
    journal_path, metadata = args
    
    try:
        # Parse journal
        data = parse_journal_file(journal_path)
        if not data:
            return False
        
        nodes = data.get('nodes', [])
        if not nodes:
            log(f"⚠️  No nodes found in {journal_path}")
            return False
        
        # Create output directory structure
        output_base = Path("aide_journals")
        run_dir = output_base / metadata['run_id']
        run_dir.mkdir(parents=True, exist_ok=True)
        
        # Create filename that includes domain and task
        safe_task = metadata['task_id'].replace('/', '_').replace('\\', '_')
        filename = f"{metadata['domain']}_{safe_task}.txt"
        output_file = run_dir / filename
        
        # Generate text content
        text_steps = []
        
        for i, node in enumerate(nodes):
            step_num = i + 1
            text_steps.append(format_step_text(node, step_num))
        
        # Create header
        header_lines = [
            "=" * 100,
            "AIDE AGENT JOURNAL ANALYSIS",
            "=" * 100,
            f"Run ID:       {metadata['run_id']}",
            f"Agent:        aide", 
            f"Domain:       {metadata['domain']}",
            f"Task:         {metadata['task_id']}",
            f"Task UUID:    {metadata['task_with_uuid']}",
            f"Logs Dir:     {metadata['logs_dir']}",
            f"Generated:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Total Steps:  {len(nodes)}",
            f"Journal Path: {journal_path}",
            "=" * 100,
            "",
            ""
        ]
        
        # Write the complete file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(header_lines))
            f.write('\n\n'.join(text_steps))
        
        log(f"✅ Processed: {metadata['run_id']} -> {metadata['domain']}/{metadata['task_id']} -> {filename}")
        return True
        
    except Exception as e:
        log(f"❌ Failed to process {journal_path}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Parse AIDE journal files from failed runs into readable text format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Parse all AIDE journals from default failed_runs directory:
  python scripts/analyze_failures/parse_aide_journals.py
  
  # Parse from custom directory:
  python scripts/analyze_failures/parse_aide_journals.py --failed-runs-dir /path/to/failed_runs
  
  # Use more concurrent processing:
  python scripts/analyze_failures/parse_aide_journals.py --concurrent 10
  
  # Dry run to see what would be processed:
  python scripts/analyze_failures/parse_aide_journals.py --dry-run

Output structure:
  aide_journals/
  ├── 2025-08-21T02-06-57-GMT_run-group_aide/
  │   ├── kaggle_uw-madison-gi-tract-image-segmentation.txt
  │   └── manual_open-problems-predict-modality.txt
  └── 2025-08-20T17-09-05-GMT_run-group_aide/
      └── kaggle_uw-madison-gi-tract-image-segmentation.txt

Each text file contains:
  - Full metadata header
  - All steps with plans, code, outputs, and analysis
  - Error information and execution details
  - Clean formatting for easy human reading
        """
    )
    
    parser.add_argument(
        "--failed-runs-dir",
        default="./failed_runs",
        help="Directory containing failed runs (default: ./failed_runs)"
    )
    parser.add_argument(
        "--output-dir",
        default="./aide_journals",
        help="Output directory for parsed journals (default: ./aide_journals)"
    )
    parser.add_argument(
        "--concurrent",
        type=int,
        default=5,
        help="Number of concurrent processing threads (default: 5)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be processed without actually doing it"
    )
    
    args = parser.parse_args()
    
    # Find AIDE journal files
    journal_files = find_aide_journal_files(args.failed_runs_dir)
    if not journal_files:
        log("❌ No AIDE journal files found")
        return 1
    
    # Dry run
    if args.dry_run:
        log("Dry run - would process these journals:")
        for i, (journal_path, metadata) in enumerate(journal_files, 1):
            safe_task = metadata['task_id'].replace('/', '_').replace('\\', '_')
            filename = f"{metadata['domain']}_{safe_task}.txt"
            print(f"  {i:2d}. {metadata['run_id']} -> {filename}")
        log(f"Total: {len(journal_files)} journal files")
        return 0
    
    # Process journals
    log(f"Processing {len(journal_files)} AIDE journal files")
    log(f"Output directory: {Path(args.output_dir).absolute()}")
    
    successful = 0
    failed = 0
    
    with ThreadPoolExecutor(max_workers=args.concurrent) as executor:
        # Submit all processing tasks
        future_to_journal = {
            executor.submit(process_journal, journal_args): journal_args 
            for journal_args in journal_files
        }
        
        # Process results
        for future in as_completed(future_to_journal):
            journal_args = future_to_journal[future]
            _, metadata = journal_args
            
            try:
                success = future.result()
                if success:
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                failed += 1
                log(f"💥 Exception processing {metadata['run_id']}: {e}")
    
    # Summary
    total = successful + failed
    log("=" * 50)
    log("PROCESSING COMPLETE")
    log(f"Total journals: {total}")
    log(f"Successful: {successful}")
    log(f"Failed: {failed}")
    log(f"Success rate: {successful/total*100:.1f}%" if total > 0 else "N/A")
    
    if successful > 0:
        log(f"Parsed journals saved to: {Path(args.output_dir).absolute()}")
        log("Text files are ready for human review!")
    
    return 0 if failed == 0 else 1

if __name__ == "__main__":
    exit(main()) 