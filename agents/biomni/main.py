"""
Runs Biomni.
"""

import argparse
import getpass
import os
import sys
from pathlib import Path

from biomni.agent import A1


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run Biomni agent")
    parser.add_argument(
        "--llm", type=str, default="claude-sonnet-4-20250514", help="LLM model to use"
    )
    parser.add_argument("--path", type=str, default=None, help="Path to agent directory")
    parser.add_argument("--timeout_seconds", type=int, default=600, help="Timeout in seconds")

    # Parse arguments
    args = parser.parse_args()

    username = getpass.getuser()

    # Check if the current user ID is 0 (root user ID on Unix-like systems)
    if os.getuid() == 0:
        print(f"You are running this script as root. Your username is '{username}'.")
    else:
        print(f"You do not have root access. Your username is {username}.")

    print("The script is being run with the following python interpreter:")
    print(sys.executable)

    # Use provided path or default to current directory
    cwd = Path(__file__).parent
    agent_path = args.path if args.path else str(cwd)

    print(f"Using LLM: {args.llm}")
    print(f"Using path: {agent_path}")
    print(f"Timeout: {args.timeout_seconds} seconds")

    # Initialize the agent with data path
    agent = A1(path=agent_path, llm=args.llm)

    prompt = """
    Build a machine learning model to solve this biomedical task. Focus on understanding the dataset structure, 
    implementing appropriate data preprocessing, selecting suitable algorithms for the task type, 
    and optimizing performance. Use appropriate evaluation metrics for the task type. Full instructions 
    can be found in the instructions.txt file.
    """

    # Run the agent
    agent.go(prompt)


if __name__ == "__main__":
    main()
