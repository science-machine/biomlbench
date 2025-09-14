#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

echo "🐳 Pulling BioML-bench Prebuilt Images"
echo "======================================"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo -e "${RED}❌ Docker is not installed or not in PATH${NC}"
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    echo -e "${RED}❌ Docker daemon is not running${NC}"
    exit 1
fi

# Mapping of prebuilt images to local agent names
declare -A IMAGE_MAP=(
    ["millerh1/biomlbench-env:v0.1a"]="biomlbench-env"
    ["millerh1/aide:v0.1a"]="aide"
    ["millerh1/biomni:v0.1a"]="biomni"
    ["millerh1/mlagentbench:v0.1a"]="mlagentbench"
    ["millerh1/stella:v0.1a"]="stella"
    ["millerh1/dummy:v0.1a"]="dummy"
)

echo -e "${YELLOW}🔽 Pulling and tagging prebuilt images...${NC}"
echo "This will download several GB of data and may take some time."
echo ""

# Track success/failure
declare -a PULLED_IMAGES=()
declare -a FAILED_IMAGES=()

for remote_image in "${!IMAGE_MAP[@]}"; do
    local_image="${IMAGE_MAP[$remote_image]}"
    
    echo -e "${BLUE}Pulling $remote_image → $local_image...${NC}"
    if docker pull "$remote_image"; then
        # Tag with local name
        if docker tag "$remote_image" "$local_image"; then
            echo -e "${GREEN}✅ Successfully pulled and tagged as $local_image${NC}"
            PULLED_IMAGES+=("$local_image")
        else
            echo -e "${RED}❌ Failed to tag $remote_image as $local_image${NC}"
            FAILED_IMAGES+=("$local_image")
        fi
    else
        echo -e "${RED}❌ Failed to pull $remote_image${NC}"
        FAILED_IMAGES+=("$local_image")
    fi
    echo ""
done

# Summary
echo -e "${YELLOW}📊 Pull Summary${NC}"
echo "==============="

if [ ${#PULLED_IMAGES[@]} -gt 0 ]; then
    echo -e "${GREEN}Successfully pulled and tagged (${#PULLED_IMAGES[@]}):"
    for image in "${PULLED_IMAGES[@]}"; do
        echo "  ✅ $image"
    done
    echo ""
fi

if [ ${#FAILED_IMAGES[@]} -gt 0 ]; then
    echo -e "${RED}Failed to pull (${#FAILED_IMAGES[@]}):"
    for image in "${FAILED_IMAGES[@]}"; do
        echo "  ❌ $image"
    done
    echo ""
    echo -e "${YELLOW}Note: You can build missing images locally using ./scripts/build_agent.sh <agent>${NC}"
else
    echo -e "${GREEN}🎉 All prebuilt images pulled and tagged successfully!${NC}"
fi

echo ""
echo "Next steps:"
echo "1. Prepare task datasets: biomlbench prepare -t polarishub/tdcommons-caco2-wang"
echo "2. Run agents: biomlbench run-agent --agent aide --task-id polarishub/tdcommons-caco2-wang"
echo ""
echo "Available agents: aide, biomni, mlagentbench, stella, dummy"
