# Invokes phusion/holy-build-box-64:2.0.1 docker container to build portable pandora.
# Inspired by Páll Melsted blog (https://pmelsted.wordpress.com/2015/10/14/building-binaries-for-bioinformatics/),
# on how he et al managed to make kallisto (Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter,
# Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519)
# portable in different linux distributions.
set -eu
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PANDORA_DIR="$(dirname "${SCRIPT_DIR}")"
PORTABLE_EXECUTABLE_BUILD_DIR="${PANDORA_DIR}/build_portable_executable"

cd $PANDORA_DIR

if [ -d "${PORTABLE_EXECUTABLE_BUILD_DIR}" ]; then
  echo "Please remove ${PORTABLE_EXECUTABLE_BUILD_DIR} before proceeding."
  exit 1
fi

sudo docker run -t -i --rm \
  -v ${PANDORA_DIR}:/io \
  phusion/holy-build-box-64:2.0.1 \
  bash /io/portable_binary_builder/build_portable_binary_core.sh
