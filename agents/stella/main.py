"""
Simple BioMLBench interface for STELLA
"""
import argparse
import os
import sys
from pathlib import Path

# Add STELLA to Python path
sys.path.insert(0, "/home/agent/ref/STELLA")

# Import stella_core module early to initialize all global components
import stella_core
from stella_core import *


def main():
    parser = argparse.ArgumentParser(description="STELLA agent for BioMLBench")
    parser.add_argument("--model", type=str, default="claude-sonnet-4")
    # FIX: Use action="store_true" for booleans, not type=bool
    parser.add_argument("--use_templates", action="store_true")
    parser.add_argument("--use_mem0", action="store_true") 
    parser.add_argument("--max_tools", type=int, default=30)
    parser.add_argument("--timeout_seconds", type=int, default=21600)
    
    args = parser.parse_args()
    
    print(f"Initializing STELLA with model: {args.model}")
    
    # Verify all core components exist after stella_core import
    print(f"🤖 Dev agent: {type(stella_core.dev_agent).__name__ if stella_core.dev_agent else 'Missing'}")
    print(f"🎯 Critic agent: {type(stella_core.critic_agent).__name__ if stella_core.critic_agent else 'Missing'}")
    print(f"🛠️ Tool creation agent: {type(stella_core.tool_creation_agent).__name__ if stella_core.tool_creation_agent else 'Missing'}")
    print(f"🧠 Models available: grok={stella_core.grok_model is not None}, gemini={stella_core.gemini_model is not None}")
    
    # Initialize STELLA (creates manager_agent and memory system)
    success = stella_core.initialize_stella(
        use_template=args.use_templates,
        use_mem0=args.use_mem0
    )
    
    if not success:
        raise RuntimeError("Failed to initialize STELLA")
    
    # Access manager_agent after initialization
    manager_agent = stella_core.manager_agent
    
    if manager_agent is None:
        raise RuntimeError("Manager agent was not created properly")
    
    print(f"✅ Manager agent: {type(manager_agent).__name__} with {len(manager_agent.tools)} tools")
    
    # Read task description
    task_description_path = "/home/data/description.md"
    with open(task_description_path, 'r') as f:
        task_description = f.read()
    
    # STELLA prompt for biomedical task execution
    prompt = f"""
Build a machine learning model to solve this biomedical task. Focus on understanding the dataset structure, 
implementing appropriate data preprocessing, selecting suitable algorithms for the task type, 
and optimizing performance. Use appropriate evaluation metrics for the task type.

Task description and data are available in /home/data/. 
Save your final predictions to /home/submission/submission.csv (or appropriate format based on sample submission).

{task_description}
"""
    
    print("Running STELLA on biomedical task...")
    result = manager_agent.run(prompt)
    print("STELLA execution completed.")
    
    return result


if __name__ == "__main__":
    main() 