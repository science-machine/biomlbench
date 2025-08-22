# AWS Deployment for BioML-bench

Automated deployment system for running BioML-bench experiments on Amazon Web Services (AWS) EC2.

## AWS vs GCP Equivalents

| GCP Component | AWS Equivalent | Notes |
|--------------|----------------|-------|
| GCP VM Image | AMI (Amazon Machine Image) | Pre-built image with biomlbench installed |
| gcloud CLI | AWS CLI + boto3 | Python SDK for AWS services |
| gcloud compute ssh | AWS Systems Manager (SSM) | Secure remote command execution |
| n2-standard-16 | m5.4xlarge | 16 vCPUs, 64GB RAM |
| g2-standard-16 (L4 GPU) | g4dn.4xlarge (T4 GPU) | Similar GPU performance tier |
| Zone (us-central1-a) | Availability Zone (us-east-1a) | Regional compute resources |

## Prerequisites

### 1. AWS CLI Setup
```bash
# Install AWS CLI
pip install awscli boto3

# Configure credentials
aws configure
# Enter: Access Key ID, Secret Access Key, Region (us-east-1), Output format (json)
```

### 2. Create Required AWS Resources

#### a. Create AMI (Amazon Machine Image)
First, create an EC2 instance with the biomlbench environment:
```bash
# Launch a base instance (Ubuntu 22.04)
aws ec2 run-instances \
  --image-id ami-0c02fb55956c7d316 \
  --instance-type m5.4xlarge \
  --key-name your-key \
  --security-group-ids sg-xxxxxx \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":500,"VolumeType":"gp3"}}]'

# SSH into instance and set up biomlbench
ssh ubuntu@<instance-ip>

# Install biomlbench and dependencies
# ... (your setup commands)

# Create AMI from the configured instance
aws ec2 create-image \
  --instance-id i-xxxxx \
  --name "biomlbench-ami-v1" \
  --description "BioML-bench environment with all dependencies"
```

#### b. Create EC2 Key Pair
```bash
aws ec2 create-key-pair \
  --key-name biomlbench-key \
  --query 'KeyMaterial' \
  --output text > biomlbench-key.pem

chmod 400 biomlbench-key.pem
```

#### c. Create Security Group
```bash
# Create security group
aws ec2 create-security-group \
  --group-name biomlbench-sg \
  --description "Security group for BioML-bench instances"

# Add rules (SSH and SSM)
aws ec2 authorize-security-group-ingress \
  --group-id sg-xxxxx \
  --protocol tcp \
  --port 22 \
  --cidr 0.0.0.0/0

# For SSM, outbound HTTPS is needed (usually allowed by default)
```

#### d. Create IAM Role and Instance Profile
```bash
# Create the setup script
cat > setup-aws-resources.sh << 'EOF'
#!/bin/bash

# Create IAM role for EC2 instances
aws iam create-role --role-name biomlbench-role --assume-role-policy-document '{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "ec2.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}'

# Attach policies for S3 and SSM access
aws iam attach-role-policy \
  --role-name biomlbench-role \
  --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess

aws iam attach-role-policy \
  --role-name biomlbench-role \
  --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore

# Create instance profile
aws iam create-instance-profile --instance-profile-name biomlbench-instance-profile

# Add role to instance profile
aws iam add-role-to-instance-profile \
  --instance-profile-name biomlbench-instance-profile \
  --role-name biomlbench-role
EOF

chmod +x setup-aws-resources.sh
./setup-aws-resources.sh
```

## Usage from Project Root

The deployment system works identically to the GCP version:

```bash
# Test the deployment (dry run)
python scripts/aws-deploy/deploy-aws.py --jobs aws-test.txt --ami ami-xxxxx --dry-run

# Run with 5 concurrent instances
python scripts/aws-deploy/deploy-aws.py --jobs aws-jobs.txt --concurrent 5 --ami ami-xxxxx

# Run CPU jobs
python scripts/aws-deploy/deploy-aws.py \
  --jobs scripts/aws-deploy/production-cpu-jobs.txt \
  --concurrent 15 \
  --machine-type cpu \
  --ami ami-xxxxx

# Run GPU jobs
python scripts/aws-deploy/deploy-aws.py \
  --jobs scripts/aws-deploy/production-gpu-jobs.txt \
  --concurrent 8 \
  --machine-type gpu \
  --ami ami-xxxxx

# Use different region
python scripts/aws-deploy/deploy-aws.py \
  --jobs aws-jobs.txt \
  --concurrent 5 \
  --ami ami-xxxxx \
  --region us-west-2
```

## Job Files

Job files are identical to GCP format - just use the same files:

```bash
# Copy existing job files
cp scripts/gcp-deploy/production-cpu-jobs.txt scripts/aws-deploy/
cp scripts/gcp-deploy/production-gpu-jobs.txt scripts/aws-deploy/
```

Format remains `agent,task_id`:
```
aide,polarishub/tdcommons-caco2-wang
biomni,proteingym-dms/A0A1I9GEU1_NEIME_Kennouche_2019
```

## What It Does

For each job, the script:

1. **Creates EC2 Instance** - m5.4xlarge (CPU) or g4dn.4xlarge (GPU) with biomlbench AMI
2. **Waits for SSM** - Ensures AWS Systems Manager agent is ready
3. **Runs Pipeline via SSM**:
   - `biomlbench run-agent --agent {agent} --task-id {task_id}`
   - `biomlbench grade --submission {run}/submission.jsonl --output-dir results/`
   - Verifies S3 uploads exist
4. **Cleans Up** - Terminates instance regardless of success/failure

## Key Differences from GCP

1. **Authentication**: Uses AWS credentials instead of gcloud auth
2. **Remote Execution**: AWS SSM instead of gcloud compute ssh
3. **Instance Naming**: EC2 tags (max 255 chars) vs GCP labels (max 63 chars)
4. **Machine Types**: Different naming but equivalent specs
5. **Retry Logic**: Same infinite retry for quota handling

## Monitoring

```bash
# Watch active instances
watch 'aws ec2 describe-instances --filters "Name=tag:Project,Values=biomlbench" "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].[InstanceId,Tags[?Key=='\''Name'\''].Value|[0],State.Name]" --output table'

# Check SSM command status
aws ssm list-commands --filter key=Status,value=InProgress

# Check S3 artifacts (same as GCP)
aws s3 ls s3://biomlbench/v1/artifacts/runs/ --recursive
aws s3 ls s3://biomlbench/v1/artifacts/grades/ --recursive
```

## Cost Comparison

| Component | GCP | AWS | Notes |
|-----------|-----|-----|-------|
| CPU Instance | n2-standard-16: ~$0.77/hr | m5.4xlarge: ~$0.77/hr | Nearly identical |
| GPU Instance | g2-standard-16 (L4): ~$1.50/hr | g4dn.4xlarge (T4): ~$1.20/hr | AWS slightly cheaper |
| Storage | 500GB SSD: ~$0.05/hr | 500GB gp3: ~$0.04/hr | Similar costs |
| **Cost per Job** | ~$0.40-$0.80 | ~$0.35-$0.70 | AWS ~10-15% cheaper |

**Examples:**
- 10 jobs with 15 concurrent instances: ~$3.50-7 total
- 100 jobs with 15 concurrent instances: ~$35-70 total
- **Max throughput**: Depends on your EC2 limits (request increases as needed)

## Troubleshooting

### SSM Connection Issues
```bash
# Check if SSM agent is installed on AMI
aws ssm describe-instance-information

# Ensure instance has internet access for SSM
# Check security group allows outbound HTTPS (443)
```

### Instance Launch Failures
```bash
# Check your service quotas
aws service-quotas get-service-quota \
  --service-code ec2 \
  --quota-code L-1216C47A  # Running On-Demand instances

# Request quota increase if needed
aws service-quotas request-service-quota-increase \
  --service-code ec2 \
  --quota-code L-1216C47A \
  --desired-value 100
```

### Missing Resources
```bash
# List all resources to verify setup
aws ec2 describe-images --owners self  # Your AMIs
aws ec2 describe-key-pairs  # Key pairs
aws ec2 describe-security-groups  # Security groups
aws iam list-instance-profiles  # Instance profiles
```

## Example Full Workflow

```bash
# 1. Set up AWS resources (one-time)
./setup-aws-resources.sh

# 2. Create your AMI (one-time)
# ... follow AMI creation steps above ...

# 3. Copy job files
cp scripts/gcp-deploy/production-*.txt scripts/aws-deploy/

# 4. Test without executing
python scripts/aws-deploy/deploy-aws.py \
  --jobs scripts/aws-deploy/production-cpu-jobs.txt \
  --ami ami-xxxxx \
  --dry-run

# 5. Run CPU jobs
python scripts/aws-deploy/deploy-aws.py \
  --jobs scripts/aws-deploy/production-cpu-jobs.txt \
  --concurrent 15 \
  --machine-type cpu \
  --ami ami-xxxxx

# 6. Monitor progress
watch 'aws ec2 describe-instances --filters "Name=tag:Project,Values=biomlbench" --output table'

# 7. Check results in S3
aws s3 ls s3://biomlbench/v1/artifacts/runs/ --recursive | tail -20
```

## Security Considerations

1. **Instance Profile**: Grants S3 access only, following least privilege
2. **SSM**: More secure than SSH - no open ports, IAM-based access
3. **Temporary Resources**: All instances auto-terminate after job completion
4. **VPC**: Consider using a dedicated VPC for isolation (optional) 