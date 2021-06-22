#!/bin/bash
# run from project root: scripts/create_archives.sh <PANDORA_VERSION>
# based on https://github.com/nzanepro/git-archive-submodules/blob/master/bin/git-archive-submodules.sh
set -eu

########################################################################################################################
# argument parsing
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters."
    echo "Usage: $0 <PANDORA_VERSION>"
    echo "Example: $0 0.8.0"
    exit 1
fi
PANDORA_VERSION="$1"
########################################################################################################################

########################################################################################################################
# configs
PANDORA_URL="https://github.com/rmcolq/pandora"
ARCHIVES_DIR="./archives"
########################################################################################################################

########################################################################################################################
# main script
ARCHIVES_DIR=$(realpath "${ARCHIVES_DIR}")
if [ -d "${ARCHIVES_DIR}" ]; then
  echo "Please remove ${ARCHIVES_DIR} before proceeding."
  exit 1
fi

mkdir -p "${ARCHIVES_DIR}"
cd "${ARCHIVES_DIR}"

TARPREFIX="pandora-${PANDORA_VERSION}"
echo "Cloning ${TARPREFIX}"
git clone --recursive --depth=1 --single-branch --branch "${PANDORA_VERSION}" "${PANDORA_URL}" "${PANDORA_VERSION}"

echo "Creating tar archive..."
cd "${PANDORA_VERSION}"
git archive --prefix="${TARPREFIX}"/ -o "${ARCHIVES_DIR}/${TARPREFIX}.tar" "${PANDORA_VERSION}"
git submodule foreach --recursive \
  "git archive --prefix=${TARPREFIX}/\${displaypath}/ HEAD > ${ARCHIVES_DIR}/tmp.tar && \
  tar --concatenate --file=${ARCHIVES_DIR}/${TARPREFIX}.tar ${ARCHIVES_DIR}/tmp.tar" > /dev/null
rm "${ARCHIVES_DIR}/tmp.tar"

echo "Compressing to tar.gz..."
gzip -9 "${ARCHIVES_DIR}/${TARPREFIX}.tar"

echo "Compressing to zip..."
cd "${ARCHIVES_DIR}" && tar xzf "${TARPREFIX}.tar.gz" && \
zip -r "${TARPREFIX}.zip" "${TARPREFIX}" > /dev/null

echo "All done!"
########################################################################################################################
