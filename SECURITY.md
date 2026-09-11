# Security policy

## Supported versions

nuFor is in active development and has not published a tagged stable release,
so only the current `main` branch receives security fixes. This will be
revised once tagged releases exist.

## Reporting a vulnerability

Report security issues privately, not in public issues. Use GitHub's private
vulnerability reporting: open the repository's Security tab and choose
"Report a vulnerability". That opens a private advisory with the maintainer
instead of a public thread.

Include as much of this as you can:

- the component and file involved
- a short description of the issue
- reproduction steps where possible
- likely impact, if known

Reports are acknowledged promptly and taken to resolution.

## Scope

Anything that reads untrusted input is in scope: case files, mesh and restart
files, and the CLI and web layer, across the Rust and Fortran code alike.
Secrets, API keys, and credentials must not be committed to the repository;
they belong in environment variables or an ignored `.env` file.