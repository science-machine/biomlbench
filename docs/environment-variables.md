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

# Mem0 API key
MEM0_API_KEY=your-mem0-api-key-here

# OpenRouter API key
OPENROUTER_API_KEY=your-openrouter-api-key-here
```



### Method 2: Using Shell Environment Variables

You can also export environment variables directly in your shell:

```bash
export OPENAI_API_KEY="your-openai-api-key-here"
export ANTHROPIC_API_KEY="your-anthropic-api-key-here"
export GEMINI_API_KEY="your-gemini-api-key-here"
export MEM0_API_KEY="your-mem0-api-key-here"
export OPENROUTER_API_KEY="your-openrouter-api-key-here"
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

