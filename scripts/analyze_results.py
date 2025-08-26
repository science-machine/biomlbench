# %%

import boto3
import os
import tarfile
import pandas as pd
import glob
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# %%

v1_path = "s3://biomlbench/v1/artifacts/"
v2_path = "s3://biomlbench/v2/artifacts/"

subfolders = [
    "runs",
    "failed_runs",
    "grades"
]

download_location_v1 = "analysis_output/v1"
download_location_v2 = "analysis_output/v2"


def download_s3_object(s3_client, bucket, key, local_path):
    """Download a single S3 object to local path and extract if it's a tar.gz file."""
    try:
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        s3_client.download_file(bucket, key, local_path)
        
        # Extract tar.gz files in place
        if local_path.endswith('.tar.gz'):
            extract_dir = os.path.dirname(local_path)
            with tarfile.open(local_path, 'r:gz') as tar:
                tar.extractall(path=extract_dir)
            # Remove the tar.gz file after extraction
            os.remove(local_path)
            return f"Downloaded and extracted: {key}"
        
        return f"Downloaded: {key}"
    except Exception as e:
        return f"Failed to download/extract {key}: {str(e)}"


def download_artifacts_parallel(s3_path, local_path, max_workers=10):
    """Download artifacts from S3 in parallel."""
    # Parse S3 path
    s3_path = s3_path.replace("s3://", "")
    bucket, prefix = s3_path.split("/", 1)
    
    # Initialize S3 client
    s3_client = boto3.client('s3')
    
    # Get list of objects to download
    paginator = s3_client.get_paginator('list_objects_v2')
    objects_to_download = []
    
    for subfolder in subfolders:
        full_prefix = f"{prefix}{subfolder}/"
        for page in paginator.paginate(Bucket=bucket, Prefix=full_prefix):
            if 'Contents' in page:
                for obj in page['Contents']:
                    key = obj['Key']
                    # Create local path maintaining folder structure
                    relative_path = key.replace(prefix, "")
                    local_file_path = os.path.join(local_path, relative_path)
                    objects_to_download.append((bucket, key, local_file_path))
    
    print(f"Found {len(objects_to_download)} objects to download")
    
    # Download objects in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_key = {
            executor.submit(download_s3_object, s3_client, bucket, key, local_path): key 
            for bucket, key, local_path in objects_to_download
        }
        
        completed = 0
        for future in as_completed(future_to_key):
            result = future.result()
            completed += 1
            if completed % 10 == 0 or completed == len(objects_to_download):
                print(f"Progress: {completed}/{len(objects_to_download)} downloads completed")
    
    print(f"All downloads completed to {local_path}")



# %%
# import shutil

# # Unlink the download_location_v1 and download_location_v2 directories
# if os.path.exists(download_location_v1):
#     shutil.rmtree(download_location_v1)
# if os.path.exists(download_location_v2):
#     shutil.rmtree(download_location_v2)

# # %%

# download_artifacts_parallel(v1_path, download_location_v1)
# download_artifacts_parallel(v2_path, download_location_v2)

# %%

# Build manifest of runs

def build_runs_manifest():
    """Build a manifest of all runs from both v1 and v2 directories."""
    runs_data = []
    
    # Process both v1 and v2 directories
    for version, base_path in [("v1", download_location_v1), ("v2", download_location_v2)]:
        runs_dir = os.path.join(base_path, "runs")

        # Find all submission.jsonl files
        submission_pattern = os.path.join(runs_dir, "*", "*", "*", "submission.jsonl")
        submission_files = glob.glob(submission_pattern)
        
        print(f"Found {len(submission_files)} submission files in {version}")
        
        for submission_path in submission_files:
            try:
                # Parse the path: runs_dir/agent/task/run_id/submission.jsonl
                rel_path = os.path.relpath(submission_path, runs_dir)
                path_parts = rel_path.split(os.sep)
                
                if len(path_parts) != 4 or path_parts[3] != "submission.jsonl":
                    continue
                    
                agent = path_parts[0]
                task = path_parts[1]
                run_id = path_parts[2]
                
                # Parse task to get domain and task_id
                # Handle known multi-part domains first
                if task.startswith('proteingym-dms-'):
                    task_domain = 'proteingym-dms'
                    task_id = task[len('proteingym-dms-'):]
                elif '-' in task:
                    task_domain = task.split('-', 1)[0]
                    task_id = task.split('-', 1)[1]
                else:
                    raise ValueError(f"Task {task} does not contain a dash")
                
                # Find corresponding run.log file
                run_log_pattern = os.path.join(runs_dir, agent, task, run_id, task_domain, "*", "run.log")
                run_log_files = glob.glob(run_log_pattern)
                
                if run_log_files:
                    log_path = run_log_files[0]  # Take the first match
                else:
                    print(f"No run.log found for {submission_path}")
                    log_path = None
                
                # Create unique identifier
                unique_run_id = f"{run_id}_{agent}_{task_id}"
                
                # Get submission.csv path from log path
                if log_path:
                    submission_csv_path = log_path.replace('run.log', 'submission/submission.csv')
                else:
                    submission_csv_path = None
                
                runs_data.append({
                    'version': version,
                    'run_id': run_id,
                    'unique_run_id': unique_run_id,
                    'agent': agent,
                    'task': task,
                    'task_domain': task_domain,
                    'task_id': task_id,
                    'submission_path': submission_path,
                    'submission_csv_path': submission_csv_path,
                    'log_path': log_path
                })
                
            except Exception as e:
                print(f"Error processing {submission_path}: {e}")
                continue
    
    # Create DataFrame
    runs_df = pd.DataFrame(runs_data)
    print(f"Built manifest with {len(runs_df)} runs")
    return runs_df

# Build the manifest
runs_manifest = build_runs_manifest()

# %%

runs_manifest


# %%


# Build manifest of grades

def build_grades_manifest():
    """Build a manifest of all grading results from both v1 and v2 directories."""
    grades_data = []
    
    # Process both v1 and v2 directories
    for version, base_path in [("v1", download_location_v1), ("v2", download_location_v2)]:
        grades_dir = os.path.join(base_path, "grades")
        
        # Find all JSON files in grades directory
        json_pattern = os.path.join(grades_dir, "*", "*", "*", "*.json")
        json_files = glob.glob(json_pattern)
        
        print(f"Found {len(json_files)} grade files in {version}")
        
        for json_path in json_files:
            try:
                # Read the JSON file
                with open(json_path, 'r') as f:
                    grade_data = json.load(f)
                
                # Extract submission_path to get run_id
                submission_path = grade_data['submission_path']
                
                # Parse run_id from submission_path using regex
                # Pattern: /home/runner/biomlbench/runs/{run_id}/...
                run_id_match = re.search(r'/runs/([^/]+)/', submission_path)
                
                run_id = run_id_match.group(1)
                
                # Parse the JSON path to get agent and task info
                rel_path = os.path.relpath(json_path, grades_dir)
                path_parts = rel_path.split(os.sep)
                    
                agent = path_parts[0]
                task = path_parts[1]
                timestamp_report = path_parts[2]
                json_filename = path_parts[3]
                
                # Get task_id from the JSON data (authoritative source)
                json_task_id = grade_data.get('task_id', '')
                
                # Parse task to get domain and task_id for consistency
                if task.startswith('proteingym-dms-'):
                    task_domain = 'proteingym-dms'
                    parsed_task_id = task[len('proteingym-dms-'):]
                elif '-' in task:
                    task_domain = task.split('-', 1)[0]
                    parsed_task_id = task.split('-', 1)[1]
                else:
                    task_domain = task
                    parsed_task_id = task
                
                # Use the parsed task_id for unique_run_id to match runs manifest
                unique_run_id = f"{run_id}_{agent}_{parsed_task_id}"
                
                # Create record with all the grade data
                grade_record = {
                    'version': version,
                    'run_id': run_id,
                    'unique_run_id': unique_run_id,
                    'agent': agent,
                    'task': task,
                    'task_domain': task_domain,
                    'task_id': parsed_task_id,  # Use parsed for consistency with runs
                    'json_task_id': json_task_id,  # Keep original from JSON for reference
                    'timestamp_report': timestamp_report,
                    'grading_output_path': json_path,
                    'json_filename': json_filename
                }
                
                # Add all fields from the JSON (excluding submission_path which we already processed)
                for key, value in grade_data.items():
                    if key != 'submission_path':  # Avoid duplication
                        grade_record[key] = value
                
                grades_data.append(grade_record)
                
            except Exception as e:
                print(f"Error processing {json_path}: {e}")
                continue
    
    # Create DataFrame
    grades_df = pd.DataFrame(grades_data)
    # Split the task_id by "/" and take the last part
    grades_df['task_id'] = grades_df['task_id'].str.split('/').str[-1]
    print(f"Built grades manifest with {len(grades_df)} records")
    return grades_df

# Build the grades manifest
grades_manifest = build_grades_manifest()


# %% 

# Try to match the grades manifest with the runs manifest
combined_df = pd.merge(runs_manifest, grades_manifest, on=['unique_run_id', "run_id", 'agent', 'task', 'task_domain', 'version'], how='inner')
runs_without_grades = runs_manifest[~runs_manifest['unique_run_id'].isin(grades_manifest['unique_run_id'])]
runs_with_nan_grades = combined_df[combined_df['score'].isna()]


# %%

# Build manifest of failed runs

def build_failed_runs_manifest():
    """Build a manifest of all failed runs from both v1 and v2 directories."""
    failed_runs_data = []
    
    # Process both v1 and v2 directories
    for version, base_path in [("v1", download_location_v1), ("v2", download_location_v2)]:
        failed_runs_dir = os.path.join(base_path, "failed_runs")
        
        # Find all run.log files in failed_runs directory
        run_log_pattern = os.path.join(failed_runs_dir, "*", "*", "*", "run.log")
        run_log_files = glob.glob(run_log_pattern)
        
        print(f"Found {len(run_log_files)} failed run files in {version}")
        
        for log_path in run_log_files:
            try:
                # Parse the path: failed_runs_dir/run_id/task_domain/task_id_uuid/run.log
                rel_path = os.path.relpath(log_path, failed_runs_dir)
                path_parts = rel_path.split(os.sep)
                
                run_id = path_parts[0]
                task_domain = path_parts[1]
                task_id_with_uuid = path_parts[2]
                
                # Extract agent from run_id
                # Pattern: {timestamp}_run-group_{agent}
                agent_match = re.search(r'_run-group_(.+)$', run_id)
                if agent_match:
                    agent = agent_match.group(1)
                else:
                    agent = None
                
                # Remove UUID from task_id (everything after the last underscore followed by hex pattern)
                # Pattern: task_id_uuid where uuid is like 6077ea77-6268-409f-812b-70148eb35324
                task_id_match = re.match(r'(.+)_[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$', task_id_with_uuid)
                if task_id_match:
                    task_id = task_id_match.group(1)
                else:
                    # Fallback: just use the whole string
                    task_id = task_id_with_uuid
                
                # Reconstruct full task name
                task = f"{task_domain}-{task_id}"
                
                # Create unique identifier
                unique_run_id = f"{run_id}_{agent}_{task_id}"
                
                failed_runs_data.append({
                    'version': version,
                    'run_id': run_id,
                    'unique_run_id': unique_run_id,
                    'agent': agent,
                    'task': task,
                    'task_domain': task_domain,
                    'task_id': task_id,
                    'submission_path': None,  # Failed runs don't have submissions
                    'log_path': log_path
                })
                
            except Exception as e:
                print(f"Error processing {log_path}: {e}")
                continue
    
    # Create DataFrame
    failed_runs_df = pd.DataFrame(failed_runs_data)
    print(f"Built failed runs manifest with {len(failed_runs_df)} runs")
    return failed_runs_df

# Build the failed runs manifest
failed_runs_manifest = build_failed_runs_manifest()

# %% 

# Check if any of the failed runs are in the runs manifest
failed_runs_in_combined_df = failed_runs_manifest[failed_runs_manifest['unique_run_id'].isin(combined_df['unique_run_id'])]
# This is empty as expected

# %%

# Next step: get the list of failed runs and screen their logs to exclude certain ones
all_failed_runs = pd.concat([failed_runs_manifest, runs_without_grades, runs_with_nan_grades])
# Types of exclusions:
#   - "Reached max steps"
#   - "The grading server failed to start"
#   - "credit balance"
#   - "quota"
#   - Run completed in ([\d.]+) seconds (where the number is less than 180)

def analyze_failed_run_logs(failed_runs_df):
    """Analyze log files to determine which failed runs should be excluded for fairness."""
    exclusion_results = []
    
    # Define exclusion patterns
    exclusion_patterns = [
        ("reached_max_steps", r"Reached max steps"),
        ("grading_server_failed", r"The grading server failed to start"),
        ("credit_balance", r"credit balance"),
        ("quota", r"quota"),
    ]
    
    print(f"Analyzing {len(failed_runs_df)} failed run logs...")
    
    for idx, row in failed_runs_df.iterrows():
        log_path = row.get('log_path')
        if not log_path or not os.path.exists(log_path):
            exclusion_results.append({
                'unique_run_id': row['unique_run_id'],
                'exclusion_reason': 'no_log_file',
                'should_exclude': False,
                'log_path': log_path
            })
            continue
        
        try:
            # Read the log file
            with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                log_content = f.read()
            
            exclusion_reason = None
            should_exclude = False
            
            # Check for each exclusion pattern
            for reason, pattern in exclusion_patterns:
                if re.search(pattern, log_content, re.IGNORECASE):
                    exclusion_reason = reason
                    should_exclude = True
                    break
            
            # Check for quick completion (less than 180 seconds)
            if not should_exclude:
                completion_match = re.search(r"Run completed in ([\d.]+) seconds", log_content)
                if completion_match:
                    completion_time = float(completion_match.group(1))
                    if completion_time < 180:
                        exclusion_reason = f"quick_completion_{completion_time}s"
                        should_exclude = True
            
            # If no exclusion reason found, it's a legitimate failure
            if not exclusion_reason:
                exclusion_reason = "legitimate_failure"
                should_exclude = False
            
            exclusion_results.append({
                'unique_run_id': row['unique_run_id'],
                'exclusion_reason': exclusion_reason,
                'should_exclude': should_exclude,
                'log_path': log_path
            })
            
        except Exception as e:
            raise e
    
    # Create DataFrame with results
    exclusion_df = pd.DataFrame(exclusion_results)
    
    # Print summary
    exclusion_counts = exclusion_df['exclusion_reason'].value_counts()
    print("\nExclusion analysis summary:")
    for reason, count in exclusion_counts.items():
        print(f"  {reason}: {count}")
    
    excluded_count = exclusion_df['should_exclude'].sum()
    legitimate_failures = len(exclusion_df) - excluded_count
    print(f"\nTotal excluded: {excluded_count}")
    print(f"Legitimate failures: {legitimate_failures}")
    
    return exclusion_df

# Analyze the failed runs
failed_run_analysis = analyze_failed_run_logs(all_failed_runs)




# %%

# Now we take the runs that actually have grades and count them.
# V1
v1_runs_with_grades = combined_df[combined_df['version'] == 'v1']
v1_runs_with_grades[['agent', 'task']].value_counts()

# V2
v2_runs_with_grades = combined_df[combined_df['version'] == 'v2']
v2_runs_with_grades[['agent', 'task']].value_counts()

# Make final df
combined_df_with_grades = pd.concat([v1_runs_with_grades, v2_runs_with_grades])
combined_df_with_grades = combined_df_with_grades[
    ((combined_df_with_grades['version'] == 'v1') & (combined_df_with_grades['task_domain'] != 'manual')) |
    ((combined_df_with_grades['version'] == 'v2') & (combined_df_with_grades['task_domain'] == 'manual'))
]
combined_df_without_grades = combined_df[combined_df['score'].isna()]
combined_df_with_grades = combined_df_with_grades[combined_df_with_grades['score'].notna()]
# Filter for the ones that have at least 2
combined_df_with_grades['agent_task'] = combined_df_with_grades['agent'] + '_' + combined_df_with_grades['task']
combined_df_with_grades['agent_task'].value_counts()
combined_df_with_grades_not_dummy = combined_df_with_grades[combined_df_with_grades['agent'] != 'dummy']
too_low_count = (combined_df_with_grades_not_dummy['agent_task'].value_counts() <= 2)
too_low_count = too_low_count[too_low_count].index
too_low_count
combined_df_with_grades_not_dummy[combined_df_with_grades_not_dummy['agent_task'].isin(too_low_count)]

# All possible agent-task combinations
all_agents = combined_df_with_grades_not_dummy['agent'].unique()
all_tasks = combined_df_with_grades_not_dummy['task'].unique()
all_agent_tasks = [f"{agent}_{task}" for agent in all_agents for task in all_tasks]
# Figure out which are missing from the combined_df_with_grades
missing_agent_tasks = set(all_agent_tasks) - set(combined_df_with_grades_not_dummy['agent_task'])
# Print the missing agent-task combinations
for agent_task in missing_agent_tasks:
    print(agent_task)



# These are the only ones missing. They are missing because they have consistent failures. We will not rerun anything.
"stella_manual-open-problems-label-projection"
"mlagentbench_manual-open-problems-predict-modality"

# %% 

# Next step: look for agent-task-combos where all replicates have a 0% leaderboard score

# Find agent-task combos where all replicates have 0% leaderboard score
zero_leaderboard_combos = combined_df_with_grades_not_dummy.groupby(['agent', 'task']).agg({
    'leaderboard_percentile': ['count', 'min', 'max']
}).reset_index()

# Flatten column names
zero_leaderboard_combos.columns = ['agent', 'task', 'count', 'min_percentile', 'max_percentile']

# Filter for combos where max leaderboard percentile is 0
zero_leaderboard_combos = zero_leaderboard_combos[zero_leaderboard_combos['max_percentile'] == 0.0]

print("Agent-task combinations where all replicates have 0% leaderboard score:")
for _, row in zero_leaderboard_combos.iterrows():
    print(f"{row['agent']}_{row['task']}: {row['count']} replicates")

# Get log paths for these zero leaderboard combinations
if len(zero_leaderboard_combos) > 0:
    print("\nLog paths for zero leaderboard combinations:")
    for _, combo_row in zero_leaderboard_combos.iterrows():
        agent = combo_row['agent']
        task = combo_row['task']
        
        # Get all runs for this agent-task combination
        matching_runs = combined_df_with_grades_not_dummy[
            (combined_df_with_grades_not_dummy['agent'] == agent) & 
            (combined_df_with_grades_not_dummy['task'] == task)
        ]
        
        print(f"\n{agent}_{task}:")
        for _, run_row in matching_runs.iterrows():
            print(f"  {run_row['submission_csv_path']}")





# %%

# Next step: regrade the label projection tasks.

import subprocess

def regrade_submission(submission_csv_path, task_id):
    """Regrade a single submission and return the JSON result."""
    try:
        submission_csv_path = "scripts/" + submission_csv_path
        cmd = ['uv', 'run', 'biomlbench', 'grade-sample', submission_csv_path, task_id]
        print(f"Regrading {submission_csv_path} with task_id {task_id}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd='/home/paperspace/biomlbench')
        
        if result.returncode == 0:
            # Extract JSON from stderr (it's after "Task report:")
            lines = result.stderr.strip().split('\n')
            
            # Find the JSON block after "Task report:"
            collecting_json = False
            json_lines = []
            
            for line in lines:
                if 'Task report:' in line:
                    collecting_json = True
                    continue
                    
                if collecting_json:
                    # Remove log prefix if present
                    if line.strip().startswith('[') and '] ' in line:
                        clean_line = line.split('] ')[-1]
                    else:
                        clean_line = line.strip()
                        
                    if clean_line:
                        json_lines.append(clean_line)
            
            if json_lines:
                json_str = '\n'.join(json_lines)
                return json.loads(json_str)
        else:
            raise Exception(f"Error regrading {submission_csv_path}: {result.stderr}")
    except Exception as e:
        raise e


from tqdm import tqdm

# Regrade zero leaderboard combinations
print("Regrading zero leaderboard combinations...")
for _, combo_row in tqdm(zero_leaderboard_combos.iterrows(), total=len(zero_leaderboard_combos)):
    agent = combo_row['agent']
    task = combo_row['task']
    
    matching_runs = combined_df_with_grades_not_dummy[
        (combined_df_with_grades_not_dummy['agent'] == agent) & 
        (combined_df_with_grades_not_dummy['task'] == task)
    ]
    
    for _, run_row in tqdm(matching_runs.iterrows(), total=len(matching_runs)):
        submission_csv_path = run_row['submission_csv_path']
        if submission_csv_path and os.path.exists(submission_csv_path):
            print(f"Regrading {agent}_{task}: {submission_csv_path}")
            result = regrade_submission(submission_csv_path, run_row['json_task_id'])
            if result:
                print(f"New score: {result.get('score')}, leaderboard_percentile: {result.get('leaderboard_percentile')}")

# %%

# Regrade null results with existing submission CSV
print("\nRegrading null results with existing submission CSV...")
runs_with_nan_grades_not_dummy = combined_df_without_grades[combined_df_without_grades['agent'] != 'dummy']

# Create a copy to update with new grades
runs_with_nan_grades_updated = combined_df_without_grades.copy()

for idx, run_row in runs_with_nan_grades_not_dummy.iterrows():
    submission_csv_path = run_row['submission_csv_path']
    if submission_csv_path and os.path.exists(submission_csv_path):
        print(f"Regrading null result: {submission_csv_path}")
        result = regrade_submission(submission_csv_path, run_row['json_task_id'])
        if result:
            print(f"New score: {result.get('score')}, leaderboard_percentile: {result.get('leaderboard_percentile')}")
            
            # Update the dataframe with new grades
            if result.get('score') is not None and not pd.isna(result.get('score')):
                for key, value in result.items():
                    if key in runs_with_nan_grades_updated.columns:
                        runs_with_nan_grades_updated.loc[idx, key] = value

# %%

# Create fixed combined dataframe
# Remove dummy runs from the original combined_df

# Combine the original successful grades with the updated null grades
combined_df_fixed_grades = pd.concat([
    combined_df_with_grades, 
    runs_with_nan_grades_updated
], ignore_index=True)

# Only keep those where score is not nan
combined_df_fixed_grades = combined_df_fixed_grades[combined_df_fixed_grades['score'].notna()]

print(f"\nFixed combined dataframe has {len(combined_df_fixed_grades)} runs with valid scores")


# %% 
# Print all the run logs for failed biomni runs that were not previously excluded.

# Filter for biomni runs that should NOT be excluded (legitimate failures)
biomni_legitimate_failures = failed_run_analysis[
    (failed_run_analysis['unique_run_id'].str.contains('biomni', case=False, na=False)) &
    (failed_run_analysis['should_exclude'] == False)
]

# Exclude any that were successfully regraded (appear in final dataset)
regraded_run_ids = set(combined_df_fixed_grades['unique_run_id'].unique())
biomni_not_regraded = biomni_legitimate_failures[
    ~biomni_legitimate_failures['unique_run_id'].isin(regraded_run_ids)
]

print(f"Found {len(biomni_legitimate_failures)} total legitimate biomni failures")
print(f"Found {len(biomni_not_regraded)} legitimate biomni failures that were NOT regraded:")
print("\nBiomni Failed Run Logs (Legitimate Failures, Not Regraded):")
print("=" * 60)

for _, row in biomni_not_regraded.iterrows():
    print(f"\nRun ID: {row['unique_run_id']}")
    print(f"Exclusion Reason: {row['exclusion_reason']}")
    print(f"Log Path: {row['log_path']}")
    print("-" * 60)



# %%

# Now the time has come for plots and tables.

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Set publication-ready style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")

# Define consistent agent ordering and display names
AGENT_ORDER = ['dummy', 'mlagentbench', 'aide', 'stella', 'biomni']
AGENT_DISPLAY_NAMES = {
    'dummy': 'Dummy',
    'mlagentbench': 'MLAgentBench', 
    'aide': 'AIDE',
    'stella': 'STELLA',
    'biomni': 'Biomni'
}

def get_agent_order(agent: str) -> int:
    """Get sort order for agents."""
    try:
        return AGENT_ORDER.index(agent.lower())
    except ValueError:
        return 999

def format_agent_name(agent: str) -> str:
    """Format agent names for display."""
    return AGENT_DISPLAY_NAMES.get(agent.lower(), agent)

def calculate_sem(values):
    """Calculate standard error of the mean."""
    if len(values) == 0:
        return 0
    return np.std(values, ddof=1) / np.sqrt(len(values))

def create_performance_by_domain_plot(df, output_path="performance_by_domain.png"):
    """Create performance by domain plot similar to visualize_results.py style."""
    
    # Add task_group column for domains
    def get_task_group(task_domain):
        domain_map = {
            'kaggle': 'Biomedical Imaging',
            'proteingym-dms': 'Protein Engineering', 
            'polarishub': 'Drug Discovery',
            'manual': 'Single Cell Omics'
        }
        return domain_map.get(task_domain, task_domain)
    
    df_plot = df.copy()
    df_plot['task_group'] = df_plot['task_domain'].apply(get_task_group)
    
    # Get sorted agent list and task groups
    agents = sorted(df_plot['agent'].unique(), key=get_agent_order)
    agent_names = [format_agent_name(a) for a in agents]
    x_labels = sorted(df_plot['task_group'].unique())
    
    # Calculate task-level means for scatter points
    task_means = df_plot.groupby(['agent', 'task', 'task_group'])['leaderboard_percentile'].mean().reset_index()
    
    # Prepare data for bar plot - aggregate across tasks within each domain
    bar_data = task_means.groupby(['agent', 'task_group']).agg({
        'leaderboard_percentile': ['mean', calculate_sem]
    }).reset_index()
    bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    
    # Create figure
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    
    # Create grouped bar plot
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent', 
                ax=ax, alpha=0.8, hue_order=agents, order=x_labels,
                edgecolor='black', linewidth=2, ci=None)
    
    # Add error bars manually
    bars = ax.patches
    n_groups = len(x_labels)
    n_agents = len(agents)
    
    for agent_idx, agent in enumerate(agents):
        for group_idx, group in enumerate(x_labels):
            matching_rows = bar_data[(bar_data['agent'] == agent) & (bar_data['task_group'] == group)]
            if len(matching_rows) > 0:
                row = matching_rows.iloc[0]
                bar_idx = agent_idx * n_groups + group_idx
                if bar_idx < len(bars):
                    bar = bars[bar_idx]
                    ax.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                               yerr=row['sem'], fmt='none', c='black', capsize=5, linewidth=2)
    
    # Overlay strip plot with task-level points
    sns.stripplot(data=task_means, x='task_group', y='leaderboard_percentile', hue='agent',
                  ax=ax, dodge=True, size=5, alpha=0.7, color='black',
                  edgecolor='white', linewidth=1, 
                  hue_order=agents, order=x_labels, legend=False)
    
    # Styling
    ax.set_xlabel('Task Domain', fontsize=16)
    ax.set_ylabel('Leaderboard Percentile', fontsize=20)
    ax.set_title('Performance Distribution by Domain', fontsize=24)
    ax.tick_params(axis='x', rotation=45, labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    
    # Fix legend to show proper names
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left', fontsize=12, title_fontsize=13)
    
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Performance by domain plot saved to {output_path}")

def calculate_subdomain_statistics(df):
    """Calculate statistics for each agent/subdomain combination."""
    
    # Add task_group column for domains
    def get_task_group(task_domain):
        domain_map = {
            'kaggle': 'Biomedical Imaging',
            'proteingym-dms': 'Protein Engineering', 
            'polarishub': 'Drug Discovery',
            'manual': 'Single Cell Omics'
        }
        return domain_map.get(task_domain, task_domain)
    
    df_stats = df.copy()
    df_stats['task_group'] = df_stats['task_domain'].apply(get_task_group)
    
    results = []
    agents = sorted(df_stats['agent'].unique(), key=get_agent_order)
    
    for task_group in sorted(df_stats['task_group'].unique()):
        group_data = df_stats[df_stats['task_group'] == task_group]
        
        for agent in agents:
            agent_data = group_data[group_data['agent'] == agent]
            
            if len(agent_data) == 0:
                continue
            
            # Calculate task-level means first
            task_means = agent_data.groupby('task').agg({
                'leaderboard_percentile': 'mean',
                'above_median': 'mean',
                'any_medal': 'mean'
            }).reset_index()
            
            # Calculate statistics across tasks
            lb_mean = task_means['leaderboard_percentile'].mean()
            lb_sem = calculate_sem(task_means['leaderboard_percentile'])
            
            above_median_mean = task_means['above_median'].mean() * 100  # Convert to percentage
            above_median_sem = calculate_sem(task_means['above_median'] * 100)
            
            any_medal_mean = task_means['any_medal'].mean() * 100  # Convert to percentage
            any_medal_sem = calculate_sem(task_means['any_medal'] * 100)
            
            results.append({
                'agent': agent,
                'task_group': task_group,
                'n_tasks': len(task_means),
                'n_replicates': len(agent_data),
                'leaderboard_mean': lb_mean,
                'leaderboard_sem': lb_sem,
                'above_median_mean': above_median_mean,
                'above_median_sem': above_median_sem,
                'any_medal_mean': any_medal_mean,
                'any_medal_sem': any_medal_sem
            })
    
    return pd.DataFrame(results)

def calculate_overall_statistics(df):
    """Calculate overall statistics across all tasks."""
    
    results = []
    agents = sorted(df['agent'].unique(), key=get_agent_order)
    
    for agent in agents:
        agent_data = df[df['agent'] == agent]
        
        # Calculate task-level means first
        task_means = agent_data.groupby('task').agg({
            'leaderboard_percentile': 'mean',
            'above_median': 'mean',
            'any_medal': 'mean'
        }).reset_index()
        
        # Calculate statistics across all tasks
        lb_mean = task_means['leaderboard_percentile'].mean()
        lb_sem = calculate_sem(task_means['leaderboard_percentile'])
        
        above_median_mean = task_means['above_median'].mean() * 100
        above_median_sem = calculate_sem(task_means['above_median'] * 100)
        
        any_medal_mean = task_means['any_medal'].mean() * 100
        any_medal_sem = calculate_sem(task_means['any_medal'] * 100)
        
        results.append({
            'agent': agent,
            'n_tasks': len(task_means),
            'n_replicates': len(agent_data),
            'leaderboard_mean': lb_mean,
            'leaderboard_sem': lb_sem,
            'above_median_mean': above_median_mean,
            'above_median_sem': above_median_sem,
            'any_medal_mean': any_medal_mean,
            'any_medal_sem': any_medal_sem
        })
    
    return pd.DataFrame(results)

def format_value_with_sem(mean, sem, is_percentage=False):
    """Format value as mean ± SEM."""
    if is_percentage:
        return f"{mean:.1f} ± {sem:.1f}"
    else:
        return f"{mean:.2f} ± {sem:.2f}"

def bold_if_best(value, is_best):
    """Wrap value in LaTeX bold if it's the best."""
    if is_best:
        return f"\\textbf{{{value}}}"
    return value

def generate_subdomain_latex_table(stats_df):
    """Generate LaTeX table for subdomain results."""
    
    latex = []
    latex.append("\\begin{table}[htbp]")
    latex.append("\\centering")
    latex.append("\\caption{BioML-bench Results by Domain}")
    latex.append("\\label{tab:bioml_results_by_domain}")
    latex.append("\\resizebox{\\textwidth}{!}{%")
    latex.append("\\begin{tabular}{llcccc}")
    latex.append("\\toprule")
    latex.append("Domain & Agent & Leaderboard Percentile & Above Median (\\%) & Any Medal (\\%) & Tasks/Replicates \\\\")
    latex.append("\\midrule")
    
    # Process each task group
    for i, task_group in enumerate(sorted(stats_df['task_group'].unique())):
        if i > 0:
            latex.append("\\midrule")
        
        group_data = stats_df[stats_df['task_group'] == task_group]
        
        # Find best values for bolding
        best_lb = group_data['leaderboard_mean'].max()
        best_above_median = group_data['above_median_mean'].max()
        best_medal = group_data['any_medal_mean'].max()
        
        for j, (_, row) in enumerate(group_data.iterrows()):
            # Format values
            lb_str = format_value_with_sem(row['leaderboard_mean'], row['leaderboard_sem'])
            above_median_str = format_value_with_sem(row['above_median_mean'], row['above_median_sem'], is_percentage=True)
            any_medal_str = format_value_with_sem(row['any_medal_mean'], row['any_medal_sem'], is_percentage=True)
            
            # Bold if best
            lb_str = bold_if_best(lb_str, abs(row['leaderboard_mean'] - best_lb) < 0.01)
            above_median_str = bold_if_best(above_median_str, abs(row['above_median_mean'] - best_above_median) < 0.01)
            any_medal_str = bold_if_best(any_medal_str, abs(row['any_medal_mean'] - best_medal) < 0.01)
            
            # Domain name only on first row
            domain_str = task_group if j == 0 else ""
            
            # Format agent name
            agent_display = format_agent_name(row['agent'])
            
            # Tasks/replicates info
            tasks_repl_str = f"{row['n_tasks']}/{row['n_replicates']}"
            
            latex.append(f"{domain_str} & {agent_display} & {lb_str} & {above_median_str} & {any_medal_str} & {tasks_repl_str} \\\\")
    
    latex.append("\\bottomrule")
    latex.append("\\end{tabular}%")
    latex.append("}")
    latex.append("\\end{table}")
    
    return '\n'.join(latex)

def generate_overall_latex_table(stats_df):
    """Generate LaTeX table for overall results."""
    
    latex = []
    latex.append("\\begin{table}[htbp]")
    latex.append("\\centering")
    latex.append("\\caption{BioML-bench Overall Results}")
    latex.append("\\label{tab:bioml_results_overall}")
    latex.append("\\resizebox{\\textwidth}{!}{%")
    latex.append("\\begin{tabular}{lcccc}")
    latex.append("\\toprule")
    latex.append("Agent & Leaderboard Percentile & Above Median (\\%) & Any Medal (\\%) & Tasks/Replicates \\\\")
    latex.append("\\midrule")
    
    # Find best values for bolding
    best_lb = stats_df['leaderboard_mean'].max()
    best_above_median = stats_df['above_median_mean'].max()
    best_medal = stats_df['any_medal_mean'].max()
    
    for _, row in stats_df.iterrows():
        # Format values
        lb_str = format_value_with_sem(row['leaderboard_mean'], row['leaderboard_sem'])
        above_median_str = format_value_with_sem(row['above_median_mean'], row['above_median_sem'], is_percentage=True)
        any_medal_str = format_value_with_sem(row['any_medal_mean'], row['any_medal_sem'], is_percentage=True)
        
        # Bold if best
        lb_str = bold_if_best(lb_str, abs(row['leaderboard_mean'] - best_lb) < 0.01)
        above_median_str = bold_if_best(above_median_str, abs(row['above_median_mean'] - best_above_median) < 0.01)
        any_medal_str = bold_if_best(any_medal_str, abs(row['any_medal_mean'] - best_medal) < 0.01)
        
        # Format agent name
        agent_display = format_agent_name(row['agent'])
        
        # Tasks/replicates info
        tasks_repl_str = f"{row['n_tasks']}/{row['n_replicates']}"
        
        latex.append(f"{agent_display} & {lb_str} & {above_median_str} & {any_medal_str} & {tasks_repl_str} \\\\")
    
    latex.append("\\bottomrule")
    latex.append("\\end{tabular}%")
    latex.append("}")
    latex.append("\\end{table}")
    
    return '\n'.join(latex)

def normalize_replicates_to_four(df, target_replicates=4):
    """
    Normalize all task-agent combinations to exactly 4 replicates.
    - If < 4 replicates: pad with zeros  
    - If > 4 replicates: randomly downsample to 4
    - If exactly 4: keep as is
    """
    
    print(f"Normalizing all task-agent combinations to {target_replicates} replicates...")
    
    normalized_data = []
    
    # Group by agent and task
    grouped = df.groupby(['agent', 'task'])
    
    for (agent, task), group in grouped:
        current_count = len(group)
        
        if current_count < target_replicates:
            # Pad with zeros
            print(f"  {agent} x {task}: {current_count} replicates → padding with {target_replicates - current_count} zeros")
            
            # Keep existing data
            normalized_data.append(group)
            
            # Add zero entries for missing replicates
            zero_template = group.iloc[0].copy()  # Use first row as template
            zero_columns = ['score', 'leaderboard_percentile', 'above_median', 'any_medal', 
                           'gold_medal', 'silver_medal', 'bronze_medal']
            
            for col in zero_columns:
                if col in zero_template.index:
                    zero_template[col] = 0.0 if col == 'score' or col == 'leaderboard_percentile' else False
            
            # Add the required number of zero rows
            for i in range(target_replicates - current_count):
                zero_row = zero_template.copy()
                # Give each zero row a unique identifier
                zero_row['unique_run_id'] = f"{zero_row['unique_run_id']}_zero_{i+1}"
                normalized_data.append(pd.DataFrame([zero_row]))
                
        elif current_count > target_replicates:
            # Randomly downsample
            print(f"  {agent} x {task}: {current_count} replicates → downsampling to {target_replicates}")
            sampled = group.sample(n=target_replicates, random_state=42)  # Fixed seed for reproducibility
            normalized_data.append(sampled)
            
        else:
            # Exactly target_replicates - keep as is
            print(f"  {agent} x {task}: {current_count} replicates → keeping all")
            normalized_data.append(group)
    
    # Concatenate all data
    normalized_df = pd.concat(normalized_data, ignore_index=True)
    
    # Verify the normalization worked
    final_counts = normalized_df.groupby(['agent', 'task']).size()
    print(f"\nNormalization complete:")
    print(f"  Total rows: {len(df)} → {len(normalized_df)}")
    print(f"  All task-agent combos now have exactly {target_replicates} replicates: {(final_counts == target_replicates).all()}")
    
    return normalized_df

# Check replicate counts before normalization
print("📊 Analyzing replicate counts before normalization...")
replicate_counts = combined_df_fixed_grades.groupby(['agent', 'task']).size().reset_index(name='replicate_count')

# Find agent-task combos with fewer than 4 replicates
fewer_than_4 = replicate_counts[replicate_counts['replicate_count'] < 4].sort_values(['agent', 'task'])

print(f"\nAgent-task combinations with fewer than 4 replicates:")
print("=" * 60)
if len(fewer_than_4) > 0:
    for _, row in fewer_than_4.iterrows():
        print(f"{row['agent']} x {row['task']}: {row['replicate_count']} replicates")
    
    print(f"\nSummary:")
    print(f"  - Total combinations: {len(replicate_counts)}")
    print(f"  - With < 4 replicates: {len(fewer_than_4)}")
    print(f"  - Will be zero-padded: {len(fewer_than_4)}")
    
    # Show breakdown by agent
    agent_breakdown = fewer_than_4.groupby('agent').size()
    print(f"\nBreakdown by agent:")
    for agent, count in agent_breakdown.items():
        print(f"  - {agent}: {count} combinations need padding")
else:
    print("✅ All agent-task combinations already have 4+ replicates!")

# Normalize replicates before analysis
print("\n🔄 Normalizing replicates...")
combined_df_normalized = normalize_replicates_to_four(combined_df_fixed_grades)

# Create output directory
output_dir = "analysis_outputs"
os.makedirs(output_dir, exist_ok=True)

# Generate the performance by domain plot
print("Creating performance by domain plot...")
create_performance_by_domain_plot(combined_df_normalized, output_path=f"{output_dir}/performance_by_domain.png")

# Calculate statistics 
print("Calculating statistics...")
subdomain_stats = calculate_subdomain_statistics(combined_df_normalized)
overall_stats = calculate_overall_statistics(combined_df_normalized)

# Generate LaTeX tables
print("Generating LaTeX tables...")
subdomain_table = generate_subdomain_latex_table(subdomain_stats)
overall_table = generate_overall_latex_table(overall_stats)

# Save tables
print("Saving tables...")
with open(f"{output_dir}/table_by_domain.tex", 'w') as f:
    f.write(subdomain_table)

with open(f"{output_dir}/table_overall.tex", 'w') as f:
    f.write(overall_table)

# Save CSV versions for inspection
subdomain_stats.to_csv(f"{output_dir}/subdomain_statistics.csv", index=False)
overall_stats.to_csv(f"{output_dir}/overall_statistics.csv", index=False)

print("✅ Analysis complete!")
print("📊 Generated files in analysis_outputs/:")
print("  - performance_by_domain.png: Performance plot by domain")
print("  - table_by_domain.tex: LaTeX table with results by domain")
print("  - table_overall.tex: LaTeX table with overall results") 
print("  - subdomain_statistics.csv: Raw statistics by domain")
print("  - overall_statistics.csv: Raw overall statistics")

print("\n📈 Final Dataset Summary:")
final_counts = combined_df_normalized.groupby(['agent', 'task']).size()
print(f"  - Total runs: {len(combined_df_normalized)}")
print(f"  - Agent-task combinations: {len(final_counts)}")
print(f"  - All combinations have exactly 4 replicates: {(final_counts == 4).all()}")
print(f"  - Zero-padded entries: {len(combined_df_normalized[combined_df_normalized['unique_run_id'].str.contains('_zero_', na=False)])}") 





