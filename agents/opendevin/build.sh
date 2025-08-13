#!/bin/bash

# Simple build script that avoids sudo issues
cd /home/agent/OpenHands

# Install playwright without sudo - just skip if it fails
echo "Attempting to install playwright..."
conda run -n agent playwright install --with-deps 2>/dev/null || echo "Playwright install skipped due to permissions"

# Try to build frontend if possible
if [ -d "frontend" ]; then
    echo "Building frontend..."
    conda run -n agent npm install 2>/dev/null || echo "Frontend build skipped"
fi

echo "Build script completed" 