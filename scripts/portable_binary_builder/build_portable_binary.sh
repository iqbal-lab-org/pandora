#!/usr/bin/env bash

set -eu
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTS_DIR="$(dirname "${CURRENT_DIR}")"
PANDORA_DIR="$(dirname "${SCRIPTS_DIR}")"
PORTABLE_EXECUTABLE_BUILD_DIR="${PANDORA_DIR}/build_portable_executable"

cd "$PANDORA_DIR"

if [ -d "${PORTABLE_EXECUTABLE_BUILD_DIR}" ]; then
  echo "Please remove ${PORTABLE_EXECUTABLE_BUILD_DIR} before proceeding."
  exit 1
fi

docker run -t -i --rm \
  -v "${PANDORA_DIR}":/pandora \
  leandroishilima/pandora_static_binary_toolchain:0.0.1 \
  bash /pandora/scripts/portable_binary_builder/build_portable_binary_core.sh
