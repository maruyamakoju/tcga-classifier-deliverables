# Security Policy

## Supported Versions

Only the latest tagged release is supported for security and integrity fixes.

## Reporting Issues

Open a private GitHub security advisory if the repository is hosted on GitHub,
or contact the maintainer directly using the Git commit author email configured
for this repository.

Do not include private patient data, controlled-access data, or unpublished
clinical information in public issues.

## Scope

Relevant issues include:

- unsafe deserialization or code execution paths
- release bundle or checksum integrity failures
- dependency changes that break the lightweight runtime contract
- misleading behavior that could cause accidental clinical or cross-platform use

This project is not intended for clinical diagnosis or patient management.
