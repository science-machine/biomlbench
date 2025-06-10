#!/bin/bash
# BioML-bench Environment Validation Script
# Tests all critical dependencies and potential failure points

set -e

echo "🧪 BioML-bench Environment Validation"
echo "====================================="

# Test 1: Python Environment
echo "1️⃣ Testing Python environment..."
python -c "import sys; print(f'Python {sys.version}')"

# Test 2: Biomedical Dependencies
echo "2️⃣ Testing biomedical dependencies..."
python -c "
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import pandas as pd
import numpy as np
import sklearn
print('✅ RDKit version:', rdkit.__version__)
print('✅ Pandas version:', pd.__version__)
print('✅ Scikit-learn version:', sklearn.__version__)
"

# Test 3: RDKit Molecular Processing
echo "3️⃣ Testing RDKit molecular processing..."
python -c "
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
# Test SMILES processing
mol = Chem.MolFromSmiles('CCO')
assert mol is not None, 'SMILES parsing failed'
# Test fingerprint generation
fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
assert len(fp) == 2048, 'Fingerprint generation failed'
# Test descriptor calculation
mw = Descriptors.MolWt(mol)
assert mw > 0, 'Descriptor calculation failed'
print('✅ RDKit molecular processing working')
"

# Test 4: File System Permissions
echo "4️⃣ Testing file system permissions..."
TEST_DIR="/tmp/biomlbench_test_$$"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# Test data file creation
echo "id,value" > test_data.csv
echo "0,1.23" >> test_data.csv

# Test submission file creation
python -c "
import pandas as pd
df = pd.DataFrame({'id': [0, 1], 'prediction': [1.23, 4.56]})
df.to_csv('test_submission.csv', index=False)
print('✅ File creation/writing working')
"

# Test file operations that agents commonly do
python -c "
import pandas as pd
import os
# Test reading
df = pd.read_csv('test_data.csv')
assert len(df) == 1, 'File reading failed'
# Test writing with different patterns
df_out = pd.DataFrame({'id': [0], 'pred': [2.34]})
df_out.to_csv('./submission.csv', index=False)
df_out.to_csv('submission.csv', index=False)
# Test directory operations
os.makedirs('./output', exist_ok=True)
df_out.to_csv('./output/submission.csv', index=False)
print('✅ Common agent file operations working')
"

# Cleanup
cd /
rm -rf "$TEST_DIR"

# Test 5: Memory and Performance
echo "5️⃣ Testing memory and performance..."
python -c "
import numpy as np
import time
# Test large array operations (simulating molecular datasets)
start = time.time()
X = np.random.random((1000, 2048))  # Simulate 1000 molecules with ECFP
y = np.random.random(1000)
# Simulate training
from sklearn.ensemble import RandomForestRegressor
model = RandomForestRegressor(n_estimators=10, random_state=42)
model.fit(X, y)
pred = model.predict(X[:10])
elapsed = time.time() - start
print(f'✅ ML operations completed in {elapsed:.2f}s')
"

# Test 6: Container Environment (if in Docker)
echo "6️⃣ Testing container environment..."
if [ -f /.dockerenv ]; then
    echo "📦 Running in Docker container"
    echo "Available disk space:"
    df -h /tmp
    echo "Memory info:"
    free -h | head -2
    echo "CPU info:"
    nproc
else
    echo "🖥️ Running on host system"
fi

echo ""
echo "🎉 All environment validation tests passed!"
echo "System is ready for BioML-bench agent execution." 