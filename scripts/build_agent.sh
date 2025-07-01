#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to show usage
show_usage() {
    echo "Usage: $0 <agent-id>"
    echo ""
    echo "Available agents:"
    for agent_dir in agents/*/; do
        if [[ -f "$agent_dir/config.yaml" ]]; then
            agent_name=$(basename "$agent_dir")
            echo "  - $agent_name"
        fi
    done
    echo ""
    echo "Examples:"
    echo "  $0 dummy          # Build the dummy agent"
    echo "  $0 aide           # Build the AIDE agent"
}

# Check arguments
if [[ $# -eq 0 ]]; then
    echo -e "${RED}❌ No agent specified${NC}"
    show_usage
    exit 1
fi

AGENT_ID="$1"

echo -e "${BLUE}🤖 Building BioML-bench Agent: $AGENT_ID${NC}"
echo "======================================="

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

# Check if biomlbench-env base image exists
if ! docker images biomlbench-env | grep -q biomlbench-env; then
    echo -e "${RED}❌ Base image 'biomlbench-env' not found${NC}"
    echo "Please build the base environment first:"
    echo "  ./scripts/build_base_env.sh"
    exit 1
fi

# Check if agent directory exists
AGENT_DIR="agents/$AGENT_ID"
if [[ ! -d "$AGENT_DIR" ]]; then
    echo -e "${RED}❌ Agent directory not found: $AGENT_DIR${NC}"
    show_usage
    exit 1
fi

echo -e "${YELLOW}📋 Pre-build checks for agent '$AGENT_ID'...${NC}"

# Check if required agent files exist
required_files=(
    "$AGENT_DIR/Dockerfile"
    "$AGENT_DIR/config.yaml"
    "$AGENT_DIR/start.sh"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo -e "${RED}❌ Required file missing: $file${NC}"
        exit 1
    fi
    echo -e "${GREEN}✅ $file${NC}"
done

# Check if the agent is properly configured
echo "Validating agent configuration..."
if uv run python -c "from agents.registry import registry; agent = registry.get_agent('$AGENT_ID'); print(f'Agent {agent.id} loaded successfully')"; then
    echo -e "${GREEN}✅ Agent configuration is valid${NC}"
else
    echo -e "${RED}❌ Agent configuration validation failed${NC}"
    exit 1
fi

echo -e "${YELLOW}🔨 Building agent Docker image...${NC}"

# Build the agent Docker image
BUILD_CMD="docker build --platform=linux/amd64 -t $AGENT_ID $AGENT_DIR/ \
    --build-arg SUBMISSION_DIR=/home/submission \
    --build-arg LOGS_DIR=/home/logs \
    --build-arg CODE_DIR=/home/code \
    --build-arg AGENT_DIR=/home/agent"

echo "Running: $BUILD_CMD"
echo ""

if eval $BUILD_CMD; then
    echo ""
    echo -e "${GREEN}✅ Successfully built agent image: $AGENT_ID${NC}"
else
    echo -e "${RED}❌ Failed to build agent image: $AGENT_ID${NC}"
    exit 1
fi

echo -e "${YELLOW}🧪 Testing agent image...${NC}"

# Test agent image creation
echo "Testing agent container creation..."
if docker create --name test-agent-$AGENT_ID $AGENT_ID > /dev/null 2>&1; then
    docker rm test-agent-$AGENT_ID > /dev/null 2>&1
    echo -e "${GREEN}✅ Agent container can be created${NC}"
else
    echo -e "${RED}❌ Agent container creation failed${NC}"
    exit 1
fi

# Test agent registry loading
echo "Testing agent registry integration..."
if uv run python -c "from agents.registry import registry; agent = registry.get_agent('$AGENT_ID'); print(f'Agent image: {agent.name}')"; then
    echo -e "${GREEN}✅ Agent registry integration works${NC}"
else
    echo -e "${RED}❌ Agent registry integration failed${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}🎉 Agent '$AGENT_ID' build completed successfully!${NC}"
echo ""
echo "Next steps:"
echo "1. Test agent execution:"
echo "   biomlbench run-agent --agent $AGENT_ID --task-list experiments/splits/caco2-wang.txt"
echo ""
echo "2. Grade agent results:"
echo "   biomlbench grade --submission submission.jsonl --output-dir results/"
echo ""
echo "Image details:"
docker images $AGENT_ID --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}" 