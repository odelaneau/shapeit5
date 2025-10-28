#!/usr/bin/env bash
# Configure this repo to use versioned git hooks from .githooks/
# Usage: bash scripts/setup-git-hooks.sh
set -euo pipefail

# Ensure we're in the repo root
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." &>/dev/null && pwd)"
cd "$REPO_ROOT"

# Set hooks path
if git rev-parse --git-dir >/dev/null 2>&1; then
  git config core.hooksPath .githooks
else
  echo "This does not appear to be a Git repository." >&2
  exit 1
fi

# Ensure hook is executable
if [ -f .githooks/pre-commit ]; then
  chmod +x .githooks/pre-commit
fi

echo "Configured Git hooks path to .githooks and ensured pre-commit is executable."

echo "Tip: You can change the size limit via MAX_GIT_FILE_SIZE_MB, e.g.:"
echo "     MAX_GIT_FILE_SIZE_MB=100 git commit -m 'Allow 100MB temporarily'"

# If Git LFS is available, ensure it's installed for this repo
if command -v git >/dev/null 2>&1 && git lfs version >/dev/null 2>&1; then
  git lfs install
  echo "Git LFS detected and initialized for this repository."
else
  echo "Git LFS not detected. For large binaries consider installing: https://git-lfs.com/" >&2
fi
