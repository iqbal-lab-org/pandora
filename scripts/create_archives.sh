#!/bin/bash
# run from project root: scripts/create_archives.sh
# based on https://github.com/nzanepro/git-archive-submodules/blob/master/bin/git-archive-submodules.sh
set -eu

ARCHIVES_DIR="./archives"

if [ -d "${ARCHIVES_DIR}" ]; then
  echo "Please remove ${ARCHIVES_DIR} before proceeding."
  exit 1
fi

TARMODULE=$(basename "$(git rev-parse --show-toplevel)")
TARVERSION=$(git describe --tags --abbrev=0)
TARPREFIX="${TARMODULE}-${TARVERSION}"

echo "Creating tar archive..."
mkdir -p archives
git archive --prefix="${TARPREFIX}"/ -o archives/"${TARPREFIX}".tar "${TARVERSION}"
git submodule foreach --recursive \
  "git archive --prefix=${TARPREFIX}/\${displaypath}/ HEAD > ../../archives/tmp.tar && \
  tar --concatenate --file=../../archives/${TARPREFIX}.tar ../../archives/tmp.tar" > /dev/null
rm archives/tmp.tar

echo "Compressing to tar.gz..."
gzip -9 archives/"${TARPREFIX}".tar

echo "Compressing to zip..."
cd archives && tar xzf "${TARPREFIX}".tar.gz && \
zip -r "${TARPREFIX}".zip "${TARPREFIX}" > /dev/null
