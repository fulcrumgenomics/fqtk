# Contributing to fqtk

Thank you for your interest in contributing to fqtk! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Issue Tracker](#issue-tracker)
- [Pull Requests](#pull-requests)
- [Coding Guidelines](#coding-guidelines)
- [Testing](#testing)
- [Documentation](#documentation)
- [Release Process](#release-process)
- [License](#license)

## Code of Conduct

We are committed to providing a welcoming and inclusive environment for all contributors. Please be respectful and considerate in all interactions.

## Getting Started

### Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) (stable toolchain)
- Git

### Setting Up Your Development Environment

1. Fork the repository on GitHub

2. Clone your fork:
   ```console
   git clone https://github.com/YOUR_USERNAME/fqtk.git
   cd fqtk
   ```

3. Add the upstream repository as a remote:
   ```console
   git remote add upstream https://github.com/fulcrumgenomics/fqtk.git
   ```

4. Build the project:
   ```console
   cargo build
   ```

5. Run the tests to ensure everything is working:
   ```console
   cargo test
   ```

## Development Workflow

### Before Making Changes

1. Sync your fork with upstream:
   ```console
   git fetch upstream
   git checkout main
   git merge upstream/main
   ```

2. Create a new branch for your changes:
   ```console
   git checkout -b feature/your-feature-name
   ```

### Making Changes

1. Make your changes in your feature branch

2. Ensure your code follows the project's coding standards (see [Coding Guidelines](#coding-guidelines))

3. Write or update tests as needed

4. Run the pre-commit checks before committing:
   ```console
   ./ci/check.sh
   ```

   This script runs:
   - `cargo fmt` - Code formatting
   - `cargo clippy` - Linting
   - `cargo test` - Unit tests

### Committing Changes

- Write clear, concise commit messages
- Use the imperative mood in commit messages (e.g., "Add feature" not "Added feature")
- Reference relevant issues in your commit messages when applicable (e.g., "Fix #123")

## Issue Tracker

### Reporting Bugs

Before submitting a bug report:

1. Search existing issues to see if the problem has already been reported
2. Try to reproduce the issue with the latest version

When submitting a bug report, please include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior vs. actual behavior
- Your environment (OS, Rust version, fqtk version)
- Any relevant error messages or logs
- Sample input files if applicable (or a description of the input format)

### Suggesting Features

Feature requests are welcome! When suggesting a feature:

1. Check existing issues to avoid duplicates
2. Clearly describe the use case and why the feature would be valuable
3. If possible, outline how you envision the feature working

### Issue Labels

- `bug` - Something isn't working
- `enhancement` - New feature or request
- `documentation` - Documentation improvements
- `good first issue` - Good for newcomers

## Pull Requests

### Before Submitting

1. Ensure all tests pass locally: `cargo test`
2. Run the full pre-commit checks: `./ci/check.sh`
3. Update documentation if needed
4. Add tests for new functionality

### Submitting a Pull Request

1. Push your branch to your fork:
   ```console
   git push origin feature/your-feature-name
   ```

2. Open a pull request against the `main` branch of the upstream repository

3. In your PR description:
   - Describe what changes you made and why
   - Reference any related issues (e.g., "Closes #123")
   - Note any breaking changes

### Review Process

- All PRs require review before merging
- CI checks must pass (formatting, linting, tests)
- Address reviewer feedback promptly
- Keep PRs focused and reasonably sized when possible

## Coding Guidelines

### Style

fqtk follows standard Rust conventions and uses automated tools to ensure consistency:

- **Formatting**: Use `cargo fmt` to format your code
- **Linting**: Use `cargo clippy` and resolve all warnings

Run both with:
```console
cargo fmt --all
cargo clippy --all-features --all-targets -- -D warnings
```

### Best Practices

- Write idiomatic Rust code
- Prefer clarity over cleverness
- Document public APIs with doc comments
- Handle errors appropriately using `Result` and `anyhow`
- Avoid `unwrap()` in library code; use proper error handling
- Keep functions focused and reasonably sized

### Dependencies

- Be conservative when adding new dependencies
- Prefer well-maintained, widely-used crates
- Justify new dependencies in your PR description

## Testing

### Running Tests

```console
# Run all tests
cargo test

# Run tests with output
cargo test -- --nocapture

# Run a specific test
cargo test test_name
```

### Writing Tests

- Add tests for new functionality
- Add regression tests for bug fixes
- Place unit tests in the same file as the code being tested using `#[cfg(test)]` modules
- Use descriptive test names that explain what is being tested

## Documentation

### Updating the README

The repository README contains a copy of the `fqtk demux` usage, which must be manually kept in sync with the tool.

When updating the tool's arguments or docstring, refresh the README by running the included script:

```console
bash .github/scripts/update-docs.sh
```

GitHub Actions will verify that the README remains up-to-date, and commits with outdated usage will fail CI.

### Code Documentation

- Add doc comments (`///`) to public functions, structs, and modules
- Include examples in doc comments where helpful
- Keep documentation up to date with code changes

## Release Process

Releases are managed by the maintainers using [`cargo-release`](https://github.com/crate-ci/cargo-release).

This project follows [Semantic Versioning](https://semver.org/):

- **MAJOR** version for incompatible API changes
- **MINOR** version for backwards-compatible new functionality
- **PATCH** version for backwards-compatible bug fixes

See the [README](README.md#releasing-a-new-version) for detailed release instructions.

## License

By contributing to fqtk, you agree that your contributions will be licensed under the [MIT License](LICENSE).

All contributions must be compatible with the MIT license. Do not include code from sources with incompatible licenses.

---

Thank you for contributing to fqtk!
