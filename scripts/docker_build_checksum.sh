#!/usr/bin/env bash

#
# Generates checksum of an amalgamation of some of the important build config files.
# If any of these files change, this checksum changes and the 'builder' container image gets updated on CI.
# The check currently runs on every push.
#

set -euo pipefail
trap "exit" INT

cat \
  .dockerignore \
  .gitignore \
  Dockerfile \
  Makefile \
  docker-dev \
  requirements.txt \
  scripts/docker_build_checksum.sh \
  setup.py \
| sha256sum | cut -f 1 -d " "
