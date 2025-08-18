# GCP Deployment for BioML-bench

Automated deployment system for running BioML-bench experiments on Google Cloud Platform.

## Usage from Project Root

The recommended way to use this deployment system is from the **project root directory**:

```bash
# Test the deployment (dry run)
python scripts/gcp-deploy/deploy.py --jobs gcp-test.txt --dry-run

# Run with 5 concurrent VMs  
python scripts/gcp-deploy/deploy.py --jobs gcp-jobs.txt --concurrent 5

# Run in different zone
python scripts/gcp-deploy/deploy.py --jobs gcp-jobs.txt --zone us-west1-b --concurrent 3
```

## Job Files

Create job files in the **project root** with `agent,task_id` format:

**`gcp-test.txt`** (minimal test):
```
dummy,polarishub/tdcommons-caco2-wang
```

**`gcp-jobs.txt`** (full job list):
```
# BioML-bench GCP Jobs
aide,polarishub/tdcommons-caco2-wang
biomni,polarishub/tdcommons-caco2-wang
aide,proteingym-dms/A0A1I9GEU1_NEIME_Kennouche_2019
dummy,manual/histopathologic-cancer-detection
```

## What It Does

For each job, the script:

1. **Creates VM** - g2-standard-8 (L4 GPU) with pre-built biomlbench image
2. **Waits for SSH** - Ensures VM is ready
3. **Runs Pipeline**:
   - `biomlbench run-agent --agent {agent} --task-id {task_id}`
   - `biomlbench grade --submission {run}/submission.jsonl --output-dir results/`
   - Verifies S3 uploads exist
4. **Cleans Up** - Deletes VM regardless of success/failure

## Key Features

- **Infinite Retry**: VM creation retries until success (handles quotas)
- **Parallel Execution**: Configurable concurrent VMs (default: 15, max: 16 based on L4 quota)
- **Auto-Cleanup**: VMs always deleted after job completion
- **S3 Verification**: Ensures artifacts uploaded before success
- **Path Resolution**: Job files resolved relative to current directory

## Prerequisites

1. **GCP CLI configured** with biomlbench project
2. **Biomlbench VM image** available in your project
3. **S3 credentials** configured on the VM image

## Monitoring

```bash
# Watch active VMs
watch 'gcloud compute instances list | grep bioml'

# Check S3 artifacts  
aws s3 ls s3://biomlbench/runs/ --recursive
aws s3 ls s3://biomlbench/grades/ --recursive
```

## Example Full Workflow

```bash
# 1. Edit your job list
vim gcp-jobs.txt

# 2. Test without executing
python scripts/gcp-deploy/deploy.py --jobs gcp-jobs.txt --dry-run

# 3. Run with 5 concurrent VMs
python scripts/gcp-deploy/deploy.py --jobs gcp-jobs.txt --concurrent 5

# 4. Monitor progress
watch 'gcloud compute instances list | grep bioml'

# 5. Check results in S3
aws s3 ls s3://biomlbench/runs/ --recursive
```

## Cost Estimation

- **VM Cost**: ~$1.50/hour per g2-standard-8
- **Disk Cost**: ~$0.05/hour per 500GB
- **Average Job**: 15-30 minutes
- **Cost per Job**: ~$0.40-$0.80

**Examples:**
- 10 jobs with 15 concurrent VMs: ~$4-8 total
- 100 jobs with 15 concurrent VMs: ~$40-80 total  
- **Max throughput**: 16 concurrent VMs = ~32-64 jobs/hour 