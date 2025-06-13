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
- **Data source integrations** - Support new biomedical benchmark databases
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
- Update documentation for user-facing changes

## Documentation

### Building Documentation

BioML-bench uses MkDocs with auto-generated API documentation. When contributing:

**Build and serve docs locally:**
```bash
# Build static docs
./scripts/build_docs.sh

# Build and serve docs locally
./scripts/serve_docs.sh
```

**Documentation is automatically generated from:**

- **Docstrings** - All API documentation comes from code docstrings
- **CLI help** - CLI documentation is extracted from argparse definitions
- **Manual content** - User guides, tutorials, and conceptual documentation

**When to update documentation:**

- Add docstrings to new functions/classes
- Update docstrings when changing function signatures
- Add examples to docstrings for complex functionality
- Manual docs are rarely needed for API changes

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