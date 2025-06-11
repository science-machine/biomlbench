# Contributing to BioML-bench

We welcome contributions to BioML-bench! This guide covers how to contribute effectively.

## Development Setup

1. **Fork and clone the repository**
2. **Install development dependencies:**
   ```bash
   uv sync --extra dev
   ```
3. **Build the base environment:**
   ```bash
   ./scripts/build_base_env.sh
   ```

## Types of Contributions

- **New biomedical tasks** - Add tasks from your domain expertise
- **Agent implementations** - Contribute new AI agents
- **Data source integrations** - Support new biomedical databases
- **Documentation improvements** - Help improve this documentation
- **Bug fixes and optimizations** - Code quality improvements

## Contribution Process

1. Create a feature branch
2. Make your changes
3. Test thoroughly with existing tasks
4. Submit a pull request
5. Address review feedback

## Code Standards

- Follow PEP 8 style guidelines
- Add comprehensive docstrings
- Include type hints
- Write tests for new functionality
- Update documentation for user-facing changes

## Testing

Test your changes with:
```bash
# Run basic tests
./scripts/test_environment.sh

# Test specific agent
biomlbench run-agent --agent dummy --task-id caco2-wang

# Validate new tasks
biomlbench prepare -t my-new-task
``` 