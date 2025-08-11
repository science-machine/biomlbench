#!/usr/bin/env python3
"""
Wrapper for AutoBioML that configures it to run without Docker.

This wrapper patches AutoGen's DockerCommandLineCodeExecutor to use 
LocalCommandLineCodeExecutor instead, avoiding Docker-in-Docker issues.
"""

import os
import sys
from pathlib import Path

# Set environment variables to disable Docker
os.environ["AUTOBIOML_NO_DOCKER"] = "1"

# Track if patching was successful
patching_successful = False

# Monkey-patch the autogen docker executor to use local executor instead
try:
    # Import the official LocalCommandLineCodeExecutor
    from autogen_ext.code_executors.local import LocalCommandLineCodeExecutor
    
    # Monkey-patch the Docker executor to use Local executor
    import autogen_ext.code_executors.docker as docker_module
    
    # Store the original class for reference
    _OriginalDockerExecutor = docker_module.DockerCommandLineCodeExecutor
    
    # Replace with LocalCommandLineCodeExecutor
    docker_module.DockerCommandLineCodeExecutor = LocalCommandLineCodeExecutor
    
    # Also patch it in the main autogen_ext.code_executors namespace if it exists there
    try:
        import autogen_ext.code_executors as executors
        executors.DockerCommandLineCodeExecutor = LocalCommandLineCodeExecutor
    except:
        pass
    
    print("✓ Successfully patched Docker executor to use Local executor")
    patching_successful = True
    
except ImportError as e:
    print(f"⚠ Could not import autogen_ext modules: {e}")
    # Fallback: try the old autogen import paths
    try:
        from autogen.coding import LocalCommandLineCodeExecutor
        import autogen_ext.code_executors.docker as docker_module
        docker_module.DockerCommandLineCodeExecutor = LocalCommandLineCodeExecutor
        print("✓ Successfully patched using legacy autogen imports")
        patching_successful = True
    except Exception as e2:
        print(f"✗ Failed to patch with legacy imports: {e2}")
except Exception as e:
    print(f"✗ Unexpected error during patching: {e}")

# If patching failed, we need to inform the user
if not patching_successful:
    print("\n⚠ WARNING: Could not patch Docker executor to use Local executor.")
    print("AutoBioML may attempt to use Docker-in-Docker which will likely fail.")
    print("Please ensure autogen-ext is properly installed in the environment.")
    print("\nContinuing anyway...\n")

# Now run autobioml with the provided arguments
if __name__ == "__main__":
    try:
        # Pass all arguments to autobioml
        from autobioml.cli import main
        main()
    except Exception as e:
        print(f"\n✗ Error running AutoBioML: {e}")
        if not patching_successful:
            print("\nThis may be due to the Docker executor patching failure.")
            print("Consider checking the autogen installation.")
        sys.exit(1) 