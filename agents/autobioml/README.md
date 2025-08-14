# AutoBioML Agent for BioML-bench

This is the integration of the AutoBioML framework into BioML-bench.

## Overview

AutoBioML is an autonomous agent framework that uses teams of specialized AI agents to collaboratively solve biomedical machine learning challenges. It features:

- **Team A (Planning)**: Principal Scientist, ML Expert, and Bioinformatics Expert who analyze tasks and create implementation plans
- **Team B (Implementation)**: Implementation Engineer and Data Science Critic who execute plans with iterative feedback

## Agent Variants

- `autobioml` - Default configuration using GPT-4o
- `autobioml/claude` - Uses Claude 3.5 Sonnet
- `autobioml/gpt-4.1` - Uses GPT-4.0-2024-08-06
- `autobioml/dev` - Development mode with fewer iterations

## Architecture

The integration consists of several components:

1. **adapter.py** - Converts BioML-bench task format to AutoBioML's expected format
2. **autobioml_wrapper.py** - Patches AutoBioML to use AutoGen's LocalCommandLineCodeExecutor instead of DockerCommandLineCodeExecutor, enabling code execution without Docker-in-Docker
3. **convert_output.py** - Converts AutoBioML's output back to BioML-bench submission format
4. **start.sh** - Orchestrates the entire process

## How It Works

1. The adapter reads the BioML-bench task and creates a compatible `task.yaml`
2. AutoBioML is executed with the wrapper to handle local code execution
3. The agent teams collaborate through multiple iterations:
   - Team A analyzes and plans
   - Team B implements and gets feedback
   - The cycle continues until completion
4. The output is converted to BioML-bench submission format

## Requirements

- OpenAI API key (for GPT models)
- Anthropic API key (for Claude variant)
- Perplexity API key (optional, for research capabilities)

## Building

To build the agent Docker image:

```bash
./scripts/build_agent.sh autobioml
```

## Running

To run the agent on a task:

```bash
biomlbench run-agent --agent autobioml --task-id manual/task-name
```

## Notes

- The agent runs code directly in the container using AutoGen's LocalCommandLineCodeExecutor (Docker-in-Docker is disabled)
- All output files are saved to the working directory
- The lab notebook provides a detailed record of the agent's reasoning process
- Code execution happens within the BioML-bench container environment, which provides isolation from the host system 