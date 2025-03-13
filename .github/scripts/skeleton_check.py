"""Skeleton structure checks for nuFor CI (bootstrap step 1).

Verifies the top-level files the repo contract requires and that tracked text
files are LF-only, so line-ending drift never sneaks in (project rule 6).
Runs on Linux and Windows runners. Exits nonzero on any failure.
"""

from __future__ import annotations

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
    "plans/PLAN.md",
    ".github/ISSUE_TEMPLATE/bug_report.md",
    ".github/ISSUE_TEMPLATE/feature_request.md",
    ".github/workflows/ci.yml",
]

BINARY_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".webp", ".pdf", ".zip",
                   ".gz", ".tar", ".h5", ".o", ".mod", ".so", ".dll", ".exe"}

TEXT_SUFFIXES = {".md", ".txt", ".rs", ".f90", ".f", ".py", ".toml", ".yaml",
                 ".yml", ".json", ".sh", ".bat", ".ps1", ".cmake", ".svg",
                 ".html", ".css", ".js", ".yml"}


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

    # Verify text files are LF-only.
    for path in root.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        # Build artifacts (Cargo target/, CMake build/) hold binary data with
        # arbitrary bytes; they are gitignored and never part of the text tree.
        if "target" in path.parts or "build" in path.parts:
            continue
        if path.name in ("Cargo.lock",) or path.suffix in BINARY_SUFFIXES:
            continue
        raw = path.read_bytes()
        if b"\r\n" in raw:
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