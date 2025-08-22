#!/bin/bash
#
# Helper script to create a BioML-bench AMI on AWS
# This guides you through the process of creating an EC2 instance and converting it to an AMI
#

set -e

echo "==========================================="
echo "BioML-bench AMI Creation Helper"
echo "==========================================="
echo ""
echo "This script will help you create an AMI with biomlbench installed."
echo ""

# Check prerequisites
if ! command -v aws &> /dev/null; then
    echo "❌ AWS CLI is not installed. Please install it first: pip install awscli"
    exit 1
fi

# Get Ubuntu 22.04 AMI ID for the current region
REGION=$(aws configure get region)
echo "Finding Ubuntu 22.04 AMI for region: $REGION"

# Query for the latest Ubuntu 22.04 LTS AMI
UBUNTU_AMI=$(aws ec2 describe-images \
    --owners 099720109477 \
    --filters \
        "Name=name,Values=ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-*" \
        "Name=state,Values=available" \
    --query "sort_by(Images, &CreationDate)[-1].ImageId" \
    --output text)

echo "Found Ubuntu AMI: $UBUNTU_AMI"
echo ""

# Load configuration if it exists
if [ -f "aws-deploy-config.txt" ]; then
    source aws-deploy-config.txt
    echo "Loaded configuration from aws-deploy-config.txt"
else
    echo "⚠️  No aws-deploy-config.txt found. Run setup-aws-resources.sh first!"
    exit 1
fi

echo ""
echo "Step 1: Launch EC2 Instance"
echo "============================"
echo ""
echo "Run this command to launch an instance:"
echo ""
cat << EOF
aws ec2 run-instances \\
    --image-id $UBUNTU_AMI \\
    --instance-type m5.4xlarge \\
    --key-name $AWS_DEPLOY_KEY_NAME \\
    --security-group-ids $AWS_DEPLOY_SECURITY_GROUP \\
    --iam-instance-profile Name=$AWS_DEPLOY_INSTANCE_PROFILE \\
    --block-device-mappings '[{
        "DeviceName": "/dev/sda1",
        "Ebs": {
            "VolumeSize": 1000,
            "VolumeType": "gp3",
            "DeleteOnTermination": true
        }
    }]' \\
    --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=biomlbench-ami-builder}]' \\
    --output table
EOF

echo ""
echo "Note the Instance ID from the output above!"
echo ""
read -p "Press Enter after you've launched the instance and noted the Instance ID..."

echo ""
echo "Step 2: Get Instance Details"
echo "============================"
echo ""
read -p "Enter the Instance ID: " INSTANCE_ID

# Get instance public IP
PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids $INSTANCE_ID \
    --query "Reservations[0].Instances[0].PublicIpAddress" \
    --output text)

echo "Instance Public IP: $PUBLIC_IP"
echo ""

echo "Step 3: SSH Setup Script"
echo "========================"
echo ""
echo "Once the instance is running, SSH into it and run these commands:"
echo ""
cat << 'EOF' > biomlbench-setup.sh
#!/bin/bash
# BioML-bench AMI Setup Script

set -e

# Update system
sudo apt-get update
sudo apt-get upgrade -y

# Install Python 3.11
sudo apt-get install -y software-properties-common
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get install -y python3.11 python3.11-venv python3.11-dev

# Install system dependencies
sudo apt-get install -y \
    git \
    curl \
    wget \
    build-essential \
    libssl-dev \
    libffi-dev \
    python3-pip \
    awscli

# Install Docker (for some tasks)
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Start Docker service and enable on boot
sudo systemctl start docker
sudo systemctl enable docker

# Add users to docker group
sudo usermod -aG docker ubuntu
sudo usermod -aG docker runner || true

# Ensure Docker daemon is accessible
sudo chmod 666 /var/run/docker.sock

# Install NVIDIA drivers and CUDA (for GPU instances)
# Skip this if creating CPU-only AMI
if lspci | grep -i nvidia; then
    echo "GPU detected, installing NVIDIA drivers..."
    # Add NVIDIA package repositories
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
    sudo dpkg -i cuda-keyring_1.1-1_all.deb
    sudo apt-get update
    
    # Install CUDA
    sudo apt-get install -y cuda-drivers
fi

# Create runner user (matching GCP setup)
sudo useradd -m -s /bin/bash runner || true
sudo usermod -aG sudo runner
echo "runner ALL=(ALL) NOPASSWD: ALL" | sudo tee /etc/sudoers.d/runner

# Add home directory to runner user
sudo -u runner mkdir -p /home/runner
sudo -u runner chmod 755 /home/runner

# Switch to runner user directory
cd /home/runner

# Clone biomlbench repository (replace with your repo if different)
sudo -u runner git clone https://github.com/science-machine/biomlbench.git

# Set up Python environment
cd biomlbench
sudo -u runner python3.11 -m venv .venv
sudo -u runner .venv/bin/pip install --upgrade pip
sudo -u runner .venv/bin/pip install -e .

# Install additional dependencies
sudo -u runner .venv/bin/pip install boto3

# Configure AWS CLI for S3 access (will use instance profile)
sudo -u runner aws configure set region us-east-1

# Install SSM agent (usually pre-installed on Amazon Linux, but not Ubuntu)
wget https://s3.amazonaws.com/ec2-downloads-windows/SSMAgent/latest/debian_amd64/amazon-ssm-agent.deb
sudo dpkg -i amazon-ssm-agent.deb
sudo systemctl enable amazon-ssm-agent
sudo systemctl start amazon-ssm-agent

# Create directory structure
sudo -u runner mkdir -p /home/runner/biomlbench/{runs,results,logs}

# Set permissions
sudo chown -R runner:runner /home/runner

echo "✅ BioML-bench setup complete!"
echo ""
echo "Test the installation:"
echo "sudo -u runner /home/runner/biomlbench/.venv/bin/biomlbench --help"
EOF

echo ""
echo "The setup script has been saved to: biomlbench-setup.sh"
echo ""
echo "Copy and run it on the instance:"
echo "scp -i $AWS_DEPLOY_KEY_NAME.pem biomlbench-setup.sh ubuntu@$PUBLIC_IP:~/"
echo "ssh -i $AWS_DEPLOY_KEY_NAME.pem ubuntu@$PUBLIC_IP"
echo "chmod +x biomlbench-setup.sh && ./biomlbench-setup.sh"
echo ""
read -p "Press Enter after the setup is complete on the instance..."

echo ""
echo "Step 4: Create AMI"
echo "=================="
echo ""
echo "Now create an AMI from the configured instance:"
echo ""
cat << EOF
aws ec2 create-image \
    --instance-id $INSTANCE_ID \
    --name "biomlbench-ami-$(date +%Y%m%d-%H%M%S)" \
    --description "BioML-bench AMI with all dependencies installed" \
    --output table
EOF

echo ""
echo "Note the AMI ID from the output above!"
echo ""
read -p "Enter the AMI ID: " AMI_ID

echo ""
echo "Step 5: Wait for AMI Creation"
echo "============================="
echo "Waiting for AMI to be available..."

aws ec2 wait image-available --image-ids $AMI_ID

echo "✅ AMI is now available!"
echo ""

echo "Step 6: Cleanup"
echo "==============="
read -p "Do you want to terminate the builder instance? (y/n): " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    aws ec2 terminate-instances --instance-ids $INSTANCE_ID
    echo "✅ Instance terminated"
fi

echo ""
echo "==========================================="
echo "AMI Creation Complete!"
echo "==========================================="
echo ""
echo "AMI ID: $AMI_ID"
echo ""
echo "You can now use this AMI with the deployment script:"
echo ""
echo "python scripts/aws-deploy/deploy-aws.py \\"
echo "  --jobs scripts/aws-deploy/aws-test.txt \\"
echo "  --ami $AMI_ID \\"
echo "  --concurrent 2"
echo ""

# Save AMI ID to config
echo "export AWS_DEPLOY_AMI=\"$AMI_ID\"" >> aws-deploy-config.txt
echo "✅ AMI ID saved to aws-deploy-config.txt" 