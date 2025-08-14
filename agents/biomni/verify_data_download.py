#!/usr/bin/env python3
"""
Verify that Biomni data lake files were pre-downloaded correctly.
This script checks if all expected files are present in the biomni_data directory.
"""

import os
import sys
from pathlib import Path

# Import from the installed Biomni package
from biomni.env_desc import data_lake_dict


def main():
    """Verify that all data lake files are present."""
    print("Verifying Biomni data lake files...")

    # Check the agent directory for biomni_data
    agent_dir = os.environ.get("AGENT_DIR", ".")
    data_path = os.path.join(agent_dir, "biomni_data")

    if not os.path.exists(data_path):
        print(f"❌ Data directory not found: {data_path}")
        return 1

    data_lake_dir = os.path.join(data_path, "data_lake")
    benchmark_dir = os.path.join(data_path, "benchmark")

    # Check data lake files
    expected_data_lake_files = list(data_lake_dict.keys())
    missing_files = []
    present_files = []

    print(f"Checking {len(expected_data_lake_files)} data lake files...")

    for filename in expected_data_lake_files:
        file_path = os.path.join(data_lake_dir, filename)
        if os.path.exists(file_path):
            present_files.append(filename)
        else:
            missing_files.append(filename)

    # Check benchmark directory structure
    benchmark_ok = False
    if os.path.isdir(benchmark_dir):
        # Check for the hle directory (as expected by Biomni)
        patient_gene_detection_dir = os.path.join(benchmark_dir, "hle")
        if os.path.isdir(patient_gene_detection_dir):
            benchmark_ok = True
        else:
            # If hle doesn't exist, check if there are any files in benchmark
            benchmark_files = os.listdir(benchmark_dir)
            if benchmark_files:
                print(f"✅ Benchmark directory contains {len(benchmark_files)} files/directories")
                benchmark_ok = True

    # Check for critical sgRNA files that some tools require
    sgRNA_files = ["sgRNA_KO_SP_human.txt", "sgRNA_KO_SP_mouse.txt"]
    sgRNA_ok = True
    missing_sgrna = []

    for sgRNA_file in sgRNA_files:
        sgRNA_path = os.path.join(data_lake_dir, sgRNA_file)
        if not os.path.exists(sgRNA_path):
            missing_sgrna.append(sgRNA_file)
            sgRNA_ok = False

    # Print results
    print(f"✅ Found {len(present_files)}/{len(expected_data_lake_files)} data lake files")

    if missing_files:
        print(f"❌ Missing {len(missing_files)} files:")
        for filename in missing_files[:10]:  # Show first 10 missing files
            print(f"   - {filename}")
        if len(missing_files) > 10:
            print(f"   ... and {len(missing_files) - 10} more")
        return 1

    if benchmark_ok:
        print("✅ Benchmark files present")
    else:
        print("❌ Benchmark files missing")
        return 1

    if sgRNA_ok:
        print("✅ sgRNA files present")
    else:
        print(f"⚠️  Missing sgRNA files: {missing_sgrna}")
        print("   Note: Some CRISPR tools may not work without these files")

    print("✅ All Biomni data files verified successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
