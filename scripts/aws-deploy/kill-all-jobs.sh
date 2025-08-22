#!/bin/bash

# Get all running biomlbench instances
INSTANCE_IDS=$(aws ec2 describe-instances --filters "Name=tag:Project,Values=biomlbench" "Name=instance-state-name,Values=running" --query "Reservations[*].Instances[*].InstanceId" --output text)

if [ -z "$INSTANCE_IDS" ]; then
    echo "No running biomlbench instances found"
    exit 0
fi

echo "Found instances: $INSTANCE_IDS"
echo "Killing all Docker containers on these instances..."

# Kill Docker containers on each instance in parallel
for instance_id in $INSTANCE_IDS; do
    {
        # Get instance IP
        public_ip=$(aws ec2 describe-instances --instance-ids $instance_id --query "Reservations[0].Instances[0].PublicIpAddress" --output text 2>/dev/null)
        
        if [ "$public_ip" != "None" ] && [ -n "$public_ip" ]; then
            echo "[$instance_id] Killing jobs on $public_ip..."
            ssh -o StrictHostKeyChecking=no -o ConnectTimeout=10 runner@$public_ip "cd /home/runner && bash biomlbench/scripts/killall_docker.sh" 2>/dev/null || echo "[$instance_id] Failed to connect or kill jobs"
        else
            echo "[$instance_id] No public IP available"
        fi
    } &
done

# Wait for all background jobs to complete
wait

echo "✅ Finished killing all Docker containers" 