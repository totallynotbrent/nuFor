#!/usr/bin/env python3
"""Repository hygiene gate.

Checks the two things that quietly rot a public repo:

1. private local paths must never be tracked. plans/, goals/, the
   local nuFor_README.md and the research scratch dir are working
   material, not product. docs/ and docs-site/ ARE product (the
   published site), so they are required to be tracked.
2. tracked text files must be LF-only, so line-ending drift never
   lands in a clone.

Runs on linux and windows runners. reads only git-tracked files, so
ignored work areas never enter a fresh checkout.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FORBIDDEN_TRACKED_PREFIXES = (
    "plans/",
    "goals/",
    "nuFor_README.md",
    "research/scratch/",
)

# a published repo carries only the dot-directories a reader needs.
# every other dot-directory is local working material.
ALLOWED_DOT_DIRECTORIES = (".github",)

REQUIRED_TRACKED = (
    "docs/research/results.md",
    "docs/index.md",
    "docs-site/quartz.config.yaml",
)

BINARY_SUFFIXES = {
    ".png", ".jpg", ".jpeg", ".gif", ".webp", ".pdf", ".zip",
    ".gz", ".tar", ".h5", ".o", ".mod", ".so", ".dll", ".exe",
}
BINARY_NAMES = {"Cargo.lock"}


def tracked_paths(root: Path) -> list[str]:
    """the paths git tracks, newline separated."""
    out = subprocess.run(
        ["git", "-C", str(root), "ls-files"],
        capture_output=True, check=True, text=True,
    )
    return [p for p in out.stdout.split("\n") if p]


def main() -> int:
    errors: list[str] = []
    tracked = tracked_paths(ROOT)

    for prefix in FORBIDDEN_TRACKED_PREFIXES:
        hits = [p for p in tracked if p == prefix or p.startswith(prefix)]
        if hits:
            errors.append(
                f"private path tracked: {prefix} ({len(hits)} file(s), "
                f"first: {hits[0]})"
            )

    for rel in tracked:
        head = rel.split("/", 1)[0]
        if (
            "/" in rel
            and head.startswith(".")
            and head not in ALLOWED_DOT_DIRECTORIES
        ):
            errors.append(f"local dot-directory tracked: {head}/ ({rel})")

    for rel in REQUIRED_TRACKED:
        if rel not in tracked:
            errors.append(f"required published path not tracked: {rel}")

    for rel in tracked:
        path = ROOT / rel
        if not path.is_file():
            continue
        name = path.name
        if name in BINARY_NAMES or path.suffix in BINARY_SUFFIXES:
            continue
        if b"\r\n" in path.read_bytes():
            errors.append(f"CRLF line endings in {rel}")

    if errors:
        for err in errors:
            print("FAIL", err)
        return 1

    print(f"OK: repo hygiene verified ({len(tracked)} tracked paths)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
