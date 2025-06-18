#!/bin/bash

# Local script to deploy MkDocs site to S3
# Usage: ./scripts/deploy-docs.sh

set -e

echo "🔨 Building MkDocs site..."
uv run mkdocs build --clean

echo "📤 Deploying to S3..."
aws s3 sync site/ s3://biomlbench-docs/ \
    --delete \
    --cache-control "max-age=86400" \
    --region us-west-2


echo "✅ Documentation deployed successfully!"
echo "🔗 Your site should be available at: http://static-html-pages.s3-website-us-west-2.amazonaws.com/biomlbench/" 