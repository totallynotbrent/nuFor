"""Skeleton structure checks for nuFor CI (bootstrap step 1).

Verifies the top-level files the repo contract requires and that tracked text
files are LF-only, so line-ending drift never sneaks in (project rule 6).
Runs on Linux and Windows runners. Exits nonzero on any failure.

Only tracked files are walked: gitignored work areas (plans/, docs-site/,
target/, .vscode/) never enter a fresh checkout and must not be required or
scanned.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REQUIRED = [
    "LICENSE",
    "README.md",
    "CONTRIBUTING.md",
    "CODE_OF_CONDUCT.md",
    "SECURITY.md",
    ".gitattributes",
    ".gitignore",
    "docs/README.md",
    ".github/ISSUE_TEMPLATE/bug_report.md",
    ".github/ISSUE_TEMPLATE/feature_request.md",
    ".github/workflows/ci.yml",
]

BINARY_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".webp", ".pdf", ".zip",
                   ".gz", ".tar", ".h5", ".o", ".mod", ".so", ".dll", ".exe"}
BINARY_NAMES = {"Cargo.lock"}


def tracked_files(root: Path) -> list[Path]:
    """the paths git tracks, so scans never touch ignored work areas."""
    out = subprocess.run(["git", "-C", str(root), "ls-files", "-z"],
                         capture_output=True, check=True, text=False)
    files = [root / p.decode("utf-8") for p in out.stdout.split(b"\0") if p]
    return [f for f in files if f.is_file()]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    errors: list[str] = []

    for rel in REQUIRED:
        if not (root / rel).exists():
            errors.append(f"missing required file: {rel}")

    if errors:
        for e in errors:
            print("FAIL", e)
        return 1

    # Verify tracked text files are LF-only.
    for path in tracked_files(root):
        if path.name in BINARY_NAMES or path.suffix in BINARY_SUFFIXES:
            continue
        if b"\r\n" in path.read_bytes():
            errors.append(f"CRLF line endings in {path.relative_to(root)}")

    # Require LICENSE to start with the GPL banner.
    lic = (root / "LICENSE").read_text(encoding="utf-8", errors="replace")
    if "GNU GENERAL PUBLIC LICENSE" not in lic:
        errors.append("LICENSE does not appear to be the GPL text")

    if errors:
        for e in errors:
            print("FAIL", e)
        return 1
    print("OK: skeleton structure and LF-only text verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())