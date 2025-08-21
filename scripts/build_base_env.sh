#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

echo "🧬 Building BioML-bench Base Environment"
echo "======================================="

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
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

echo -e "${YELLOW}📋 Pre-build checks...${NC}"

# Check if required files exist
required_files=(
    "environment/Dockerfile"
    "environment/requirements.txt"
    "environment/grading_server.py"
    "environment/entrypoint.sh"
    "environment/instructions.txt"
    "pyproject.toml"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo -e "${RED}❌ Required file missing: $file${NC}"
        exit 1
    fi
    echo -e "${GREEN}✅ $file${NC}"
done

echo -e "${YELLOW}🔨 Building biomlbench-env Docker image...${NC}"
echo "This may take several minutes as it installs biomedical dependencies."
echo ""

# Build the Docker image (from root directory with environment/Dockerfile)
if docker build --platform=linux/amd64 -t biomlbench-env -f environment/Dockerfile .; then
    echo ""
    echo -e "${GREEN}✅ Successfully built biomlbench-env image${NC}"
else
    echo -e "${RED}❌ Failed to build biomlbench-env image${NC}"
    exit 1
fi

echo -e "${YELLOW}🧪 Testing base image...${NC}"

# Test basic functionality (bypass entrypoint for testing)
echo "Testing Python availability..."
if docker run --rm --entrypoint python biomlbench-env --version; then
    echo -e "${GREEN}✅ Python is available${NC}"
else
    echo -e "${RED}❌ Python test failed${NC}"
    exit 1
fi

echo "Testing biomlbench import..."
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n biomlb python -c "import biomlbench; print('BioML-bench imported successfully')"; then
    echo -e "${GREEN}✅ BioML-bench is importable${NC}"
else
    echo -e "${RED}❌ BioML-bench import failed${NC}"
    exit 1
fi

echo "Testing biomedical dependencies..."
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "import rdkit; import Bio; import sklearn; print('Dependencies OK')"; then
    echo -e "${GREEN}✅ Biomedical and ML dependencies are available${NC}"
else
    echo -e "${RED}❌ Dependencies test failed${NC}"
    exit 1
fi

echo "Testing agent conda environment..."
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "print('Agent environment OK')"; then
    echo -e "${GREEN}✅ Agent conda environment is ready${NC}"
else
    echo -e "${RED}❌ Agent conda environment test failed${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}🎉 Base environment build completed successfully!${NC}"
echo ""
echo "Next steps:"
echo "1. Build agent images: ./scripts/build_agent.sh dummy"
echo "2. Test agent execution: biomlbench run-agent --agent dummy --task-list experiments/splits/caco2-wang.txt"
echo ""
echo "Image details:"
docker images biomlbench-env --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}" 