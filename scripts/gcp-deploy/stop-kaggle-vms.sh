#!/bin/bash

# Script to stop all kaggle VMs
# Usage: ./stop-kaggle-vms.sh [--force]

set -e

FORCE=false
if [[ "$1" == "--force" ]]; then
    FORCE=true
fi

echo "🔍 Finding all kaggle VMs..."

# Get all kaggle instances with their zones
INSTANCES=$(gcloud compute instances list --filter="name~.*kaggle.*" --format="value(name,zone)" 2>/dev/null)

if [[ -z "$INSTANCES" ]]; then
    echo "✅ No kaggle VMs found"
    exit 0
fi

echo ""
echo "Found the following kaggle VMs:"
echo "=================================="
while IFS=$'\t' read -r name zone; do
    if [[ -n "$name" && -n "$zone" ]]; then
        STATUS=$(gcloud compute instances describe "$name" --zone="$zone" --format="value(status)" 2>/dev/null || echo "UNKNOWN")
        echo "  $name (zone: $zone, status: $STATUS)"
    fi
done <<< "$INSTANCES"

echo ""

if [[ "$FORCE" != "true" ]]; then
    read -p "❓ Do you want to stop all these VMs? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Aborted"
        exit 0
    fi
fi

echo ""
echo "🛑 Stopping kaggle VMs in parallel..."

# Arrays to track jobs and results
declare -a PIDS
declare -a VM_NAMES
declare -a VM_ZONES

# Start all stop commands in parallel
while IFS=$'\t' read -r name zone; do
    if [[ -n "$name" && -n "$zone" ]]; then
        echo "  Starting stop for $name in zone $zone..."
        
        # Run in background and capture PID
        (
            if gcloud compute instances delete "$name" --zone="$zone" --quiet 2>/dev/null; then
                echo "STOP_SUCCESS:$name"
            else
                echo "STOP_FAILED:$name"
            fi
        ) &
        
        PIDS+=($!)
        VM_NAMES+=("$name")
        VM_ZONES+=("$zone")
    fi
done <<< "$INSTANCES"

echo "  Started ${#PIDS[@]} parallel stop operations..."
echo "  Waiting for all operations to complete..."

# Wait for all background jobs to complete and collect results
SUCCESS_COUNT=0
FAIL_COUNT=0

for i in "${!PIDS[@]}"; do
    pid=${PIDS[$i]}
    name=${VM_NAMES[$i]}
    
    # Wait for this specific PID and capture its output
    if wait $pid; then
        echo "    ✅ Successfully stopped $name"
        ((SUCCESS_COUNT++))
    else
        echo "    ❌ Failed to stop $name"
        ((FAIL_COUNT++))
    fi
done

echo ""
echo "================================================"
echo "Summary:"
echo "  Successfully stopped: $SUCCESS_COUNT VMs"
echo "  Failed to stop: $FAIL_COUNT VMs"
echo "================================================"

if [[ $FAIL_COUNT -eq 0 ]]; then
    echo "✅ All kaggle VMs stopped successfully!"
    exit 0
else
    echo "⚠️  Some VMs failed to stop. Check the output above."
    exit 1
fi 