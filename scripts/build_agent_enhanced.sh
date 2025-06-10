#!/bin/bash
# Enhanced BioML-bench Agent Build Script
# Integrates comprehensive validation and testing

set -e

AGENT_NAME="$1"
SKIP_TESTS="${2:-false}"

if [ -z "$AGENT_NAME" ]; then
    echo "Usage: $0 <agent_name> [skip_tests]"
    echo "Example: $0 aide"
    echo "Example: $0 aide true  # Skip tests"
    exit 1
fi

echo "🚀 Enhanced BioML-bench Agent Build: $AGENT_NAME"
echo "=================================================="

# Step 1: Pre-build validation
echo "1️⃣ Pre-build validation..."
if [ ! -d "agents/$AGENT_NAME" ]; then
    echo "❌ Agent directory not found: agents/$AGENT_NAME"
    exit 1
fi

required_files=("Dockerfile" "config.yaml" "start.sh")
for file in "${required_files[@]}"; do
    if [ ! -f "agents/$AGENT_NAME/$file" ]; then
        echo "❌ Required file missing: agents/$AGENT_NAME/$file"
        exit 1
    fi
done

echo "✅ Pre-build validation passed"

# Step 2: Environment validation
echo "2️⃣ Environment validation..."
if [ -f "scripts/validate_environment.sh" ]; then
    chmod +x scripts/validate_environment.sh
    if ! ./scripts/validate_environment.sh; then
        echo "❌ Environment validation failed"
        exit 1
    fi
else
    echo "⚠️ Environment validation script not found, skipping..."
fi

# Step 3: Build base environment (if needed)
echo "3️⃣ Checking base environment..."
if ! docker image inspect biomlbench-env:latest >/dev/null 2>&1; then
    echo "🔨 Base environment not found, building..."
    if ! ./scripts/build_base_env.sh; then
        echo "❌ Base environment build failed"
        exit 1
    fi
else
    echo "✅ Base environment exists"
fi

# Step 4: Build agent with enhanced error handling
echo "4️⃣ Building agent Docker image..."

# Create build context with validation
BUILD_START=$(date +%s)
BUILD_LOG="build_${AGENT_NAME}_$(date +%Y%m%d_%H%M%S).log"

if ! ./scripts/build_agent.sh "$AGENT_NAME" 2>&1 | tee "$BUILD_LOG"; then
    echo "❌ Agent build failed. Check log: $BUILD_LOG"
    exit 1
fi

BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))
echo "✅ Agent built successfully in ${BUILD_TIME}s"

# Step 5: Post-build validation
echo "5️⃣ Post-build validation..."

# Check image size and details
if command -v docker >/dev/null 2>&1; then
    IMAGE_SIZE=$(docker images --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}" | grep "^$AGENT_NAME" | awk '{print $3}' | head -1)
    echo "📊 Image size: $IMAGE_SIZE"
    
    # Warn if image is too large
    if [[ "$IMAGE_SIZE" == *"GB"* ]]; then
        SIZE_NUM=$(echo "$IMAGE_SIZE" | sed 's/GB//')
        if (( $(echo "$SIZE_NUM > 30" | bc -l) )); then
            echo "⚠️ Warning: Image is very large (${IMAGE_SIZE}). Consider optimization."
        fi
    fi
fi

# Step 6: Runtime testing (if not skipped)
if [ "$SKIP_TESTS" != "true" ]; then
    echo "6️⃣ Runtime testing..."
    
    if [ -f "scripts/test_agent_pipeline.py" ]; then
        echo "🧪 Running comprehensive agent tests..."
        TEST_START=$(date +%s)
        
        if python scripts/test_agent_pipeline.py "$AGENT_NAME"; then
            TEST_END=$(date +%s)
            TEST_TIME=$((TEST_END - TEST_START))
            echo "✅ All tests passed in ${TEST_TIME}s"
        else
            echo "❌ Agent tests failed. Agent may not work correctly."
            echo "💡 You can skip tests with: $0 $AGENT_NAME true"
            exit 1
        fi
    else
        echo "⚠️ Test framework not found, performing basic validation..."
        
        # Basic container test
        if docker run --rm "$AGENT_NAME:latest" python -c "print('Container test: OK')"; then
            echo "✅ Basic container test passed"
        else
            echo "❌ Basic container test failed"
            exit 1
        fi
    fi
else
    echo "6️⃣ Runtime testing skipped"
fi

# Step 7: Generate build report
echo "7️⃣ Generating build report..."

BUILD_REPORT="build_report_${AGENT_NAME}_$(date +%Y%m%d_%H%M%S).json"
cat > "$BUILD_REPORT" << EOF
{
  "agent_name": "$AGENT_NAME",
  "build_timestamp": "$(date -Iseconds)",
  "build_time_seconds": $BUILD_TIME,
  "image_size": "$IMAGE_SIZE",
  "tests_run": $([ "$SKIP_TESTS" != "true" ] && echo "true" || echo "false"),
  "build_success": true,
  "build_log": "$BUILD_LOG",
  "docker_version": "$(docker --version)",
  "system_info": {
    "os": "$(uname -s)",
    "arch": "$(uname -m)",
    "kernel": "$(uname -r)"
  }
}
EOF

echo "📊 Build report saved: $BUILD_REPORT"

# Step 8: Cleanup and optimization suggestions
echo "8️⃣ Cleanup and optimization..."

# Clean up old build logs (keep last 5)
ls -t build_${AGENT_NAME}_*.log 2>/dev/null | tail -n +6 | xargs rm -f 2>/dev/null || true

# Clean up old test results (keep last 5)
ls -t test_results_${AGENT_NAME}_*.json 2>/dev/null | tail -n +6 | xargs rm -f 2>/dev/null || true

echo "✅ Cleanup completed"

# Step 9: Success summary
echo ""
echo "🎉 Enhanced build completed successfully!"
echo "=================================================="
echo "Agent: $AGENT_NAME"
echo "Build time: ${BUILD_TIME}s"
echo "Image size: $IMAGE_SIZE"
echo "Build log: $BUILD_LOG"
echo "Build report: $BUILD_REPORT"
echo ""
echo "Next steps:"
echo "1. Test agent execution:"
echo "   biomlbench run-agent --agent $AGENT_NAME/dev --task-list test_taskset.txt"
echo ""
echo "2. Run full evaluation:"
echo "   biomlbench run-agent --agent $AGENT_NAME --task-list experiments/splits/caco2-wang.txt"
echo ""
echo "3. Grade results:"
echo "   biomlbench grade --submission submission.jsonl --output-dir results/"
echo ""

# Optional: Send notification (if configured)
if command -v notify-send >/dev/null 2>&1; then
    notify-send "BioML-bench Build" "Agent $AGENT_NAME built successfully in ${BUILD_TIME}s"
fi 