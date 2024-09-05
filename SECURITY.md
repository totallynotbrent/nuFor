# Security policy

## Supported versions

nuFor is in active development and has not published a stable release. As of
now only the current state of `main` receives security fixes. This policy will
be revised once tagged releases exist.

## Reporting a vulnerability

Please report security issues privately rather than in public issues. Report by
opening a GitHub issue with the `security` label if you are comfortable, or
reach out directly to the maintainers through a private channel listed on the
repository.

Please include:

- The component and file involved
- A minimal description of the issue
- Reproduction steps where possible
- Impact, if known

We aim to respond to reports promptly and will acknowledge receipt.

## Scope

Rust and Fortran runtimes, the CLI and web layer, and any code that reads
untrusted input (case files, meshes, restart files) are in scope for security
review. Do not commit secrets, API keys, or credentials to the repository.