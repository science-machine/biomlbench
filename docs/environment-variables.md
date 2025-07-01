# Environment Variables in BioML-bench

This document explains how environment variables are handled in BioML-bench and how to set them up properly.

## Overview

BioML-bench uses environment variables to securely manage API keys and configuration settings. The system has been designed to automatically load environment variables from a `.env` file in your project root directory.

## Setting Up Environment Variables

### Method 1: Using a .env File (Recommended)

1. **Create a .env file** in the root directory of your BioML-bench project:
   ```bash
   touch .env
   ```

2. **Add your API keys** to the .env file:
   ```bash
   # OpenAI API key for agents like AIDE
   OPENAI_API_KEY=your-openai-api-key-here
   
   # Anthropic API key for Claude-based agents
   ANTHROPIC_API_KEY=your-anthropic-api-key-here
   
   # Google Gemini API key
   GEMINI_API_KEY=your-gemini-api-key-here
   
   # Kaggle API credentials (optional)
   KAGGLE_USERNAME=your-kaggle-username
   KAGGLE_KEY=your-kaggle-key
   
   # Container execution settings
   I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS=false
   ```

3. **Secure your .env file** (make sure it's in .gitignore):
   ```bash
   echo ".env" >> .gitignore
   ```

### Method 2: Using Shell Environment Variables

You can also export environment variables directly in your shell:

```bash
export OPENAI_API_KEY="your-openai-api-key-here"
export ANTHROPIC_API_KEY="your-anthropic-api-key-here"
export GEMINI_API_KEY="your-gemini-api-key-here"
```

## How Environment Variables Are Loaded

BioML-bench automatically loads environment variables in the following order:

1. **CLI Startup**: When you run any `biomlbench` command, the system automatically loads the `.env` file from your project root directory.

2. **Agent Initialization**: When agents are loaded, the system also attempts to load environment variables if they haven't been loaded already.

3. **Shell Environment**: Any variables already set in your shell environment take precedence over those in the `.env` file.

## Agent Configuration

Agent configurations use a special syntax to reference environment variables:

```yaml
# In agent config files (e.g., agents/aide/config.yaml)
env_vars:
  OPENAI_API_KEY: ${{ secrets.OPENAI_API_KEY }}
  ANTHROPIC_API_KEY: ${{ secrets.ANTHROPIC_API_KEY }}
```

The `${{ secrets.VARIABLE_NAME }}` syntax tells BioML-bench to:
1. Look for the environment variable `VARIABLE_NAME`
2. Load it from the `.env` file if not already set
3. Substitute the actual value into the agent configuration

## Required Environment Variables

Different agents require different API keys:

### AIDE Agent
- `OPENAI_API_KEY`: Required for all AIDE variants

### Claude-based Agents
- `ANTHROPIC_API_KEY`: Required for agents using Claude models

### Gemini-based Agents
- `GEMINI_API_KEY`: Required for agents using Google Gemini models

### Optional Variables
- `KAGGLE_USERNAME` and `KAGGLE_KEY`: Required only for Kaggle datasets
- `I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS`: Set to "true" to allow privileged containers

## Troubleshooting

### Error: "Environment variable `OPENAI_API_KEY` is not set!"

This error occurs when:
1. You haven't created a `.env` file
2. The variable is not set in your `.env` file
3. The variable name is misspelled

**Solution:**
1. Create a `.env` file in your project root
2. Add the required API key: `OPENAI_API_KEY=your-actual-api-key`
3. Restart your biomlbench command

### Error: "Import 'dotenv' could not be resolved"

This is usually a linter warning and doesn't affect functionality. The `python-dotenv` package is included in the BioML-bench dependencies.

**Solution:**
- Ensure you're using the correct Python environment
- Run `uv sync` to install all dependencies

### .env File Not Being Loaded

If your `.env` file exists but variables aren't being loaded:

1. **Check file location**: The `.env` file must be in the project root directory (same level as `biomlbench/` folder)
2. **Check file format**: Ensure variables are in the format `VARIABLE_NAME=value` (no spaces around =)
3. **Check file permissions**: Ensure the file is readable

## Best Practices

1. **Never commit API keys**: Always add `.env` to your `.gitignore` file
2. **Use descriptive names**: Follow the pattern `SERVICE_API_KEY` for API keys
3. **Keep it simple**: Don't add spaces around the equals sign in `.env` files
4. **Test your setup**: Run a simple command like `biomlbench --help` to ensure environment loading works

## Example .env File

```bash
# BioML-bench Environment Variables
# Fill in your actual API keys below

# OpenAI API key for agents like AIDE
OPENAI_API_KEY=sk-your-openai-api-key-here

# Anthropic API key for Claude-based agents
ANTHROPIC_API_KEY=sk-ant-your-anthropic-api-key-here

# Google Gemini API key
GEMINI_API_KEY=your-gemini-api-key-here

# Kaggle API credentials (optional)
KAGGLE_USERNAME=your-kaggle-username
KAGGLE_KEY=your-kaggle-key

# Container execution settings
I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS=false
```

## Implementation Details

The environment variable loading system works as follows:

1. **CLI Entry Point** (`biomlbench/cli.py`): Automatically loads `.env` file at startup
2. **Agent Utils** (`agents/utils.py`): Provides fallback loading and helpful error messages
3. **Agent Registry** (`agents/registry.py`): Processes `${{ secrets.* }}` patterns in agent configs

This multi-layered approach ensures that environment variables are available whenever and wherever they're needed in the system. 