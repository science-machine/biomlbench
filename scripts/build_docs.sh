#!/bin/bash
# Simple script to build docs with fresh CLI documentation

set -e

echo "Generating CLI documentation..."
python scripts/extract_cli_help.py

echo "Building documentation..."
mkdocs build "$@" 