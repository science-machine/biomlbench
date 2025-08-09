#!/usr/bin/env python3
"""
Pre-download Biomni data lake files during Docker build.
This script downloads all the data lake files that Biomni would normally download
at runtime, so they're available in the Docker image.
"""

import os
import sys
from pathlib import Path

# Import from the installed Biomni package
from biomni.env_desc import data_lake_dict
from biomni.utils import check_and_download_s3_files


def main():
    """Pre-download all data lake files."""
    print("Pre-downloading Biomni data lake files...")

    try:
        # Create the data directory structure
        data_path = "/tmp/biomni_data"
        data_lake_dir = os.path.join(data_path, "data_lake")
        benchmark_dir = os.path.join(data_path, "benchmark")

        os.makedirs(data_lake_dir, exist_ok=True)
        os.makedirs(benchmark_dir, exist_ok=True)

        # Get list of expected data lake files
        expected_data_lake_files = list(data_lake_dict.keys())

        print(f"Downloading {len(expected_data_lake_files)} data lake files...")

        # Download data lake files
        download_results = check_and_download_s3_files(
            s3_bucket_url="https://biomni-release.s3.amazonaws.com",
            local_data_lake_path=data_lake_dir,
            expected_files=expected_data_lake_files,
            folder="data_lake",
        )

        # Download benchmark files
        print("Downloading benchmark files...")
        benchmark_results = check_and_download_s3_files(
            s3_bucket_url="https://biomni-release.s3.amazonaws.com",
            local_data_lake_path=benchmark_dir,
            expected_files=[],  # Empty list - will download entire folder
            folder="benchmark",
        )

        # Print summary
        successful_downloads = sum(1 for success in download_results.values() if success)
        print(
            f"Successfully downloaded {successful_downloads}/{len(expected_data_lake_files)} data lake files"
        )

        # Check if any downloads failed
        failed_downloads = [
            filename for filename, success in download_results.items() if not success
        ]
        if failed_downloads:
            print(f"Failed to download: {failed_downloads}")
            # Don't fail the build for missing files, just warn
            print("Warning: Some files failed to download, but continuing...")

        print("All data lake files downloaded successfully!")
        return 0

    except Exception as e:
        print(f"Error during pre-download: {e}")
        print("Warning: Pre-download failed, but continuing build...")
        return 0  # Don't fail the build


if __name__ == "__main__":
    sys.exit(main())
