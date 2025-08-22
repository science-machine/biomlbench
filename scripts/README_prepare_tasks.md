# BioML-bench Task Preparation Script

This script (`prepare_all_tasks.py`) automates the preparation of all BioML-bench tasks in parallel, significantly reducing the time needed to set up the complete benchmark suite.

## Features

- **Parallel Execution**: Prepare multiple tasks simultaneously using configurable thread count
- **Robust Error Handling**: Comprehensive logging and error reporting
- **Progress Tracking**: Real-time status updates and detailed completion summary
- **Flexible Input**: Support for both job files and custom task lists
- **Retry Support**: Failed tasks are saved to a file for easy retry
- **Full CLI Integration**: Supports all `biomlbench prepare` options

## Quick Start

```bash
# Prepare all tasks with default settings (4 threads)
python scripts/prepare_all_tasks.py

# Prepare with more parallelism
python scripts/prepare_all_tasks.py --threads 8

# Prepare with custom options
python scripts/prepare_all_tasks.py --threads 6 --keep-raw --skip-verification
```

## Usage Examples

### Basic Usage

```bash
# Use default settings (4 threads, production-jobs.txt)
python scripts/prepare_all_tasks.py

# Specify number of threads
python scripts/prepare_all_tasks.py --threads 8
```

### Custom Job Files

```bash
# Use a different jobs file
python scripts/prepare_all_tasks.py --jobs my-custom-jobs.txt --threads 6

# Prepare specific tasks from a list file
python scripts/prepare_all_tasks.py --task-list specific_tasks.txt --threads 4
```

### Advanced Options

```bash
# Keep raw data and overwrite existing checksums
python scripts/prepare_all_tasks.py --keep-raw --overwrite-checksums --threads 8

# Use custom data directory and skip verification
python scripts/prepare_all_tasks.py --data-dir /custom/data/path --skip-verification --threads 6

# Full customization
python scripts/prepare_all_tasks.py \
    --jobs scripts/gcp-deploy/production-jobs.txt \
    --threads 10 \
    --keep-raw \
    --overwrite-checksums \
    --overwrite-leaderboard \
    --data-dir /custom/data/path
```

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--jobs` | Path to jobs file | `scripts/gcp-deploy/production-jobs.txt` |
| `--task-list` | Path to file with specific task IDs (one per line) | None |
| `--threads` | Number of parallel threads | 4 |
| `--keep-raw` | Keep raw downloaded data | False |
| `--overwrite-checksums` | Overwrite existing checksums | False |
| `--overwrite-leaderboard` | Overwrite existing leaderboard | False |
| `--skip-verification` | Skip checksum verification | False |
| `--data-dir` | Custom data directory | Default biomlbench location |

## Input File Formats

### Jobs File Format
The script expects a jobs file with lines in the format:
```
agent,task_id
```

Example:
```
aide,polarishub/polaris-adme-fang-hppb-1
biomni,polarishub/polaris-adme-fang-hppb-1
mlagentbench,kaggle/histopathologic-cancer-detection
stella,manual/open-problems-predict-modality
```

### Task List Format
For the `--task-list` option, provide one task ID per line:
```
polarishub/polaris-adme-fang-hppb-1
kaggle/histopathologic-cancer-detection
manual/open-problems-predict-modality
```

## Performance Recommendations

### Thread Count Guidelines
- **CPU-bound tasks**: Use 1-2x the number of CPU cores
- **I/O-bound tasks**: Can use 2-4x the number of CPU cores
- **Memory constraints**: Reduce threads if tasks are memory-intensive
- **Network limitations**: Consider bandwidth when downloading large datasets

### Typical Performance
- **Small tasks** (< 100MB): 30-120 seconds each
- **Medium tasks** (100MB-1GB): 2-10 minutes each  
- **Large tasks** (> 1GB): 10-60 minutes each

### Resource Usage
- Each thread may use 500MB-2GB RAM depending on the task
- Network bandwidth is shared across all threads
- Disk I/O can be a bottleneck for large datasets

## Output and Logging

### Console Output
The script provides real-time progress updates:
```
2025-08-22 00:01:50,069 - INFO - Found 23 unique tasks to prepare
2025-08-22 00:01:50,070 - INFO - Starting parallel preparation of 23 tasks using 4 threads
2025-08-22 00:01:52,145 - INFO - Starting preparation of task: kaggle/histopathologic-cancer-detection
2025-08-22 00:03:15,234 - INFO - ✓ Successfully prepared kaggle/histopathologic-cancer-detection in 83.1s
```

### Log File
All output is also saved to `task_preparation.log` for later review.

### Summary Report
At completion, you'll see a detailed summary:
```
============================================================
PREPARATION SUMMARY
============================================================
Total tasks: 23
Successful: 21
Failed: 2

✓ Successfully prepared tasks:
  - kaggle/histopathologic-cancer-detection
  - manual/open-problems-predict-modality
  [... more tasks ...]

✗ Failed tasks:
  - polarishub/problematic-task: Connection timeout
  - manual/another-task: Missing dependencies
```

### Failed Tasks
Failed task IDs are automatically saved to `failed_tasks.txt` for easy retry:
```bash
# Retry only the failed tasks
python scripts/prepare_all_tasks.py --task-list failed_tasks.txt --threads 2
```

## Error Handling

### Common Issues
1. **Network timeouts**: Increase timeout or retry with fewer threads
2. **Memory errors**: Reduce thread count
3. **Disk space**: Ensure sufficient space for all datasets
4. **API limits**: Some tasks may have rate limits

### Recovery Strategies
1. **Partial failures**: Use the generated `failed_tasks.txt` to retry
2. **Resource constraints**: Reduce `--threads` and retry
3. **Network issues**: Add `--skip-verification` for faster retries
4. **Disk space**: Use `--keep-raw=false` to save space

## Integration with BioML-bench

This script is fully compatible with the standard `biomlbench prepare` command and uses the same:
- Configuration files
- Data directories
- Validation procedures
- Error handling

After preparation, you can run agents normally:
```bash
# All tasks are now prepared and ready
biomlbench run-agent --agent aide --task-id polarishub/polaris-adme-fang-hppb-1
```

## Best Practices

1. **Start small**: Test with a few tasks first using `--task-list`
2. **Monitor resources**: Watch CPU, memory, and disk usage
3. **Use logging**: Check `task_preparation.log` for detailed information
4. **Plan for failures**: Some tasks may fail due to external dependencies
5. **Optimize threads**: Find the sweet spot for your hardware and network
6. **Keep raw data**: Use `--keep-raw` for debugging and re-preparation

## Troubleshooting

### Script won't start
- Ensure you're in the biomlbench project root
- Verify Python environment has required packages
- Check file permissions on the script

### All tasks failing
- Verify `biomlbench prepare` works manually for one task
- Check network connectivity
- Ensure sufficient disk space

### Memory issues
- Reduce `--threads` parameter
- Monitor system memory usage
- Consider preparing tasks in smaller batches

### Network timeouts
- Reduce `--threads` to decrease concurrent downloads
- Check network stability
- Use `--skip-verification` if checksums are problematic 