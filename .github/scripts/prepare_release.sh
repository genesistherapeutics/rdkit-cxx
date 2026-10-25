#!/usr/bin/env bash
# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

# Sets Cargo.toml to the version the next release ships, commits the change if
# there is one, and writes the version and HEAD's tree to GITHUB_OUTPUT. Fails
# while Cargo.toml sets publish = false.
#
# A version with no tag yet ships as it is, so the first release of an RDKit
# series is the .0 set by the pull request that moved to it. Otherwise the patch
# number goes up by one: the first two numbers name the RDKit series, which
# build.rs checks.

set -euo pipefail

package="$(cargo metadata --no-deps --format-version 1 |
  jq -c '.packages[] | select(.name == "rdkit-cxx")')"
if [[ "$(jq -c .publish <<< "${package}")" == "[]" ]]; then
  echo "::error::Cargo.toml sets publish = false; remove it to release."
  exit 1
fi
current="$(jq -r .version <<< "${package}")"

status=0
git ls-remote --exit-code --tags origin "refs/tags/v${current}" > /dev/null || status=$?
case "${status}" in
  0)
    IFS=. read -r major minor patch <<< "${current}"
    version="${major}.${minor}.$((patch + 1))"
    # Replaces the first match only: the package's version, above any table.
    CURRENT="${current}" VERSION="${version}" perl -0pi -e \
      's/^version = "\Q$ENV{CURRENT}\E"$/version = "$ENV{VERSION}"/m' Cargo.toml
    cargo update --workspace
    git -c user.name=github-actions -c user.email=github-actions@github.com \
      commit --quiet -m "Release ${version}" Cargo.toml Cargo.lock
    ;;
  2)
    version="${current}"
    ;;
  *)
    echo "::error::git ls-remote failed with status ${status}."
    exit 1
    ;;
esac

echo "Releasing ${version}"
{
  echo "version=${version}"
  echo "tree=$(git rev-parse 'HEAD^{tree}')"
} >> "${GITHUB_OUTPUT}"
