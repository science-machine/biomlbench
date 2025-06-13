#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧪 Testing BioML-bench Environment${NC}"
echo "===================================="

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

# Check if biomlbench-env image exists
if ! docker images biomlbench-env | grep -q biomlbench-env; then
    echo -e "${RED}❌ Base image 'biomlbench-env' not found${NC}"
    echo "Please build the base environment first:"
    echo "  ./scripts/build_base_env.sh"
    exit 1
fi

echo -e "${YELLOW}🔍 Testing base environment...${NC}"

# Test 1: Basic container startup (bypass entrypoint)
echo "1. Testing container startup..."
if docker run --rm --entrypoint echo biomlbench-env "Container startup test" > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Container starts successfully${NC}"
else
    echo -e "${RED}❌ Container startup failed${NC}"
    exit 1
fi

# Test 2: Python environment
echo "2. Testing Python environment..."
python_version=$(docker run --rm --entrypoint python biomlbench-env --version 2>&1)
if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✅ Python available: $python_version${NC}"
else
    echo -e "${RED}❌ Python test failed${NC}"
    exit 1
fi

# Test 3: Conda environments
echo "3. Testing conda environments..."
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env env list | grep -q "agent"; then
    echo -e "${GREEN}✅ Agent conda environment exists${NC}"
else
    echo -e "${RED}❌ Agent conda environment missing${NC}"
    exit 1
fi

if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env env list | grep -q "biomlb"; then
    echo -e "${GREEN}✅ BioMLB conda environment exists${NC}"
else
    echo -e "${RED}❌ BioMLB conda environment missing${NC}"
    exit 1
fi

# Test 4: BioML-bench package
echo "4. Testing BioML-bench package..."
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n biomlb python -c "import biomlbench; print(f'BioML-bench version: {biomlbench.__version__}')" 2>/dev/null; then
    echo -e "${GREEN}✅ BioML-bench package is importable${NC}"
else
    echo -e "${RED}❌ BioML-bench package import failed${NC}"
    exit 1
fi

# Test 5: Biomedical dependencies
echo "5. Testing biomedical dependencies..."

# Test RDKit
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "from rdkit import Chem; print('RDKit OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ RDKit is available${NC}"
else
    echo -e "${RED}❌ RDKit test failed${NC}"
    exit 1
fi

# Test BioPython
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "import Bio; print('BioPython OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ BioPython is available${NC}"
else
    echo -e "${RED}❌ BioPython test failed${NC}"
    exit 1
fi

# Skip Polaris test since we removed it from dependencies
echo -e "${YELLOW}⏭️  Polaris test skipped (not in dependencies)${NC}"

# Test 6: Scientific computing libraries
echo "6. Testing scientific computing libraries..."

# Test numpy
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "import numpy; print('numpy OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ numpy is available${NC}"
else
    echo -e "${RED}❌ numpy test failed${NC}"
    exit 1
fi

# Test pandas
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "import pandas; print('pandas OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ pandas is available${NC}"
else
    echo -e "${RED}❌ pandas test failed${NC}"
    exit 1
fi

# Test scikit-learn (import as sklearn)
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "import sklearn; print('scikit-learn OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ scikit-learn is available${NC}"
else
    echo -e "${RED}❌ scikit-learn test failed${NC}"
    exit 1
fi

# Test scipy
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n agent python -c "import scipy; print('scipy OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ scipy is available${NC}"
else
    echo -e "${RED}❌ scipy test failed${NC}"
    exit 1
fi

# Test 7: Directory structure
echo "7. Testing container directory structure..."
expected_dirs=("/home" "/home/submission" "/private" "/opt/conda")
for dir in "${expected_dirs[@]}"; do
    if docker run --rm --entrypoint test biomlbench-env -d "$dir" 2>/dev/null; then
        echo -e "${GREEN}✅ Directory exists: $dir${NC}"
    else
        echo -e "${RED}❌ Directory missing: $dir${NC}"
        exit 1
    fi
done

# Test 8: Flask and grading server
echo "8. Testing Flask availability..."
if docker run --rm --entrypoint /opt/conda/bin/conda biomlbench-env run -n biomlb python -c "import flask; print('Flask OK')" 2>/dev/null; then
    echo -e "${GREEN}✅ Flask is available in biomlb environment${NC}"
else
    echo -e "${RED}❌ Flask test failed${NC}"
    exit 1
fi

# Test 9: Entrypoint script
echo "9. Testing entrypoint script..."
if docker run --rm --entrypoint test biomlbench-env -x /entrypoint.sh 2>/dev/null; then
    echo -e "${GREEN}✅ Entrypoint script is executable${NC}"
else
    echo -e "${RED}❌ Entrypoint script test failed${NC}"
    exit 1
fi

# Test 10: File permissions
echo "10. Testing file permissions..."
if docker run --rm --entrypoint test biomlbench-env -w /home 2>/dev/null; then
    echo -e "${GREEN}✅ /home directory is writable${NC}"
else
    echo -e "${RED}❌ /home directory permissions test failed${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}🎉 All environment tests passed!${NC}"
echo ""
echo "Environment summary:"
echo "- Base image: biomlbench-env"
echo "- Agent environment: /opt/conda/envs/agent"
echo "- Grading environment: /opt/conda/envs/biomlb"  
echo "- Biomedical libraries: RDKit, BioPython, Polaris"
echo "- Scientific computing: NumPy, Pandas, scikit-learn, SciPy"
echo ""
echo "Ready for agent builds and testing!" 