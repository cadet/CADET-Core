#!/usr/bin/env bash
set -euo pipefail

export HOMEBREW_NO_AUTO_UPDATE=1
export HOMEBREW_NO_INSTALL_CLEANUP=1

need() { brew list "$1" >/dev/null 2>&1 || brew install "$1"; }

# Build tools
need cmake
need ninja
need pkg-config

# CADET deps
need eigen@3
need hdf5
need tbb
need suite-sparse
need superlu

brew --prefix eigen
ls -la "$(brew --prefix eigen)/share/eigen3/cmake" || true
pkg-config --list-all | grep -E 'hdf5|tbb|cholmod|umfpack|superlu' || true
