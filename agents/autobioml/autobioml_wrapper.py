#!/usr/bin/env python3
"""
Wrapper script to run AutoBioML with local code execution instead of Docker.
This monkey-patches AutoGen's DockerCommandLineCodeExecutor to use LocalCommandLineCodeExecutor.
"""

import sys
import os

def patch_docker_executor():
    """Monkey-patch AutoGen's Docker executor to use Local executor instead."""
    try:
        # Import the modules we need to patch
        from autogen_ext.code_executors.local import LocalCommandLineCodeExecutor
        import autogen_ext.code_executors.docker as docker_module
        
        # Create a wrapper class that filters out Docker-specific parameters
        class FilteredLocalExecutor(LocalCommandLineCodeExecutor):
            def __init__(self, **kwargs):
                # Filter out Docker-specific parameters
                docker_specific_params = {
                    'image', 'container_name', 'auto_remove', 'stop_container', 
                    'device_requests', 'volumes', 'environment', 'network_mode',
                    'ports', 'privileged', 'working_dir', 'user'
                }
                filtered_kwargs = {k: v for k, v in kwargs.items() if k not in docker_specific_params}
                super().__init__(**filtered_kwargs)
        
        # Monkey-patch the DockerCommandLineCodeExecutor
        original_docker_executor = getattr(docker_module, 'DockerCommandLineCodeExecutor', None)
        
        if original_docker_executor:
            # Replace DockerCommandLineCodeExecutor with our filtered local executor
            docker_module.DockerCommandLineCodeExecutor = FilteredLocalExecutor
            print("✓ Successfully patched Docker executor to use Local executor")
        else:
            print("⚠ DockerCommandLineCodeExecutor not found - continuing without patch")
            
    except Exception as e:
        print(f"⚠ Failed to patch Docker executor: {e}")
        print("Continuing with original AutoBioML configuration...")

if __name__ == "__main__":
    # Apply the patch before importing AutoBioML
    patch_docker_executor()
    
    # Now import and run AutoBioML
    from autobioml.cli import autobioml
    
    # Run AutoBioML with the provided arguments
    autobioml() 