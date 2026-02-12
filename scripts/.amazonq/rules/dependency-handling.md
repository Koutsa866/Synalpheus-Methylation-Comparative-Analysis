# Dependency Handling Rule

When writing code that requires external tools or packages:
- NEVER gracefully handle missing dependencies by skipping functionality
- ALWAYS attempt to automatically install missing dependencies
- Use appropriate package managers:
  - `brew` for macOS system tools
  - `apt` for Linux system tools
  - `pip` for Python packages
- Only fail with clear error message if automatic installation is impossible
- Verify installation succeeded before proceeding
