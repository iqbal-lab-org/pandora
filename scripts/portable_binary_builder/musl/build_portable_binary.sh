#!/usr/bin/env bash

set -eu
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PANDORA_DIR="$(realpath "${CURRENT_DIR}/../../..")"
PORTABLE_EXECUTABLE_BUILD_DIR="${PANDORA_DIR}/build_portable_executable_musl"

cd "$PANDORA_DIR"

if [ -d "${PORTABLE_EXECUTABLE_BUILD_DIR}" ]; then
  echo "Please remove ${PORTABLE_EXECUTABLE_BUILD_DIR} before proceeding."
  exit 1
fi

docker run -t -i --rm \
  -v "${PANDORA_DIR}":/pandora \
  leandroishilima/pandora_static_binary_toolchain_musl:0.0.2 \
  bash /pandora/scripts/portable_binary_builder/musl/build_portable_binary_core.sh
