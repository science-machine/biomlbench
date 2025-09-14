# Cloud Deployment

BioML-bench provides automated deployment scripts for running large-scale experiments on AWS and GCP.

## Google Cloud Platform (GCP)

### Prerequisites
- GCP CLI configured: `gcloud auth login`
- Biomlbench VM image available in your project
- S3 credentials for artifact uploads

### Usage

```bash
# Create job file with agent,task_id format
cat > my-jobs.txt << EOF
dummy,polarishub/tdcommons-caco2-wang
aide,proteingym-dms/SPIKE_SARS2_Starr_2020_binding
EOF

# Run deployment (automatically creates/manages VMs)
python deploy/gcp-deploy/deploy.py \
  --jobs my-jobs.txt \
  --s3-prefix s3://your-bucket/artifacts \
  --concurrent 5
```

Features:
- **Automated**: Creates VMs, runs jobs, uploads to S3, cleans up
- **Parallel**: Configurable concurrent VMs (max 16 for L4 GPU quota)
- **Resilient**: Infinite retry for VM creation, auto-cleanup on failure

Cost: ~$0.40-0.80 per job (g2-standard-8 with L4 GPU)

## Amazon Web Services (AWS)

### Prerequisites

```bash
# Set up AWS resources (one-time)
./deploy/aws-deploy/setup-aws-resources.sh

# This creates security groups, key pairs, and saves config
source aws-deploy-config.txt
```

### Usage

```bash
# Run deployment (replace ami-xxxxx with your AMI)
python deploy/aws-deploy/deploy-aws.py \
  --jobs my-jobs.txt \
  --s3-prefix s3://your-bucket/artifacts \
  --ami ami-xxxxx \
  --security-group $AWS_DEPLOY_SECURITY_GROUP \
  --concurrent 5
```

Features:
- **Equivalent to GCP**: Same job format and workflow
- **SSM-based**: Secure remote execution without SSH
- **Cost-effective**: ~10-15% cheaper than GCP

Cost: ~$0.35-0.70 per job (m5.4xlarge CPU / g4dn.4xlarge GPU)

## Job Files

Both platforms use identical job file format:

```
# Comments allowed
agent,task_id
aide,polarishub/tdcommons-caco2-wang
biomni,proteingym-dms/SPIKE_SARS2_Starr_2020_binding
dummy,kaggle/histopathologic-cancer-detection
```

Pre-made job files:
- `deploy/gcp-deploy/production-cpu-jobs.txt` - Non-image tasks (80 jobs)
- `deploy/gcp-deploy/production-gpu-jobs.txt` - Image analysis tasks (16 jobs)

See respective deployment folder READMEs for detailed setup and usage instructions.
