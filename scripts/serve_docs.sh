#!/bin/bash
# Simple script to serve docs with fresh CLI documentation

set -e

echo "Generating CLI documentation..."
python scripts/extract_cli_help.py

echo "Serving documentation..."
mkdocs serve "$@" 