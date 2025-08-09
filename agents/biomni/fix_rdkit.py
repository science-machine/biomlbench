#!/usr/bin/env python3
"""
Fix RDKit compatibility issues at runtime.
This script downgrades NumPy if needed to ensure RDKit works.
"""

import subprocess
import sys

def fix_rdkit_compatibility():
    """Ensure RDKit is compatible with the current environment."""
    try:
        # Try importing RDKit
        from rdkit import Chem
        print("RDKit is already working correctly")
        return True
    except (ImportError, AttributeError) as e:
        print(f"RDKit import failed: {e}")
        print("Attempting to fix by downgrading NumPy...")
        
        # Downgrade NumPy
        subprocess.run([
            sys.executable, "-m", "pip", "install", "numpy<2.0", "--force-reinstall"
        ], check=True)
        
        # Try again
        try:
            from rdkit import Chem
            print("RDKit fixed successfully!")
            return True
        except Exception as e:
            print(f"Failed to fix RDKit: {e}")
            return False

if __name__ == "__main__":
    success = fix_rdkit_compatibility()
    sys.exit(0 if success else 1) 