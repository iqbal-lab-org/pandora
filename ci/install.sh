#!/bin/bash
set -evu

BOOST_URL="http://sourceforge.net/projects/boost/files/boost/${BOOST_VERSION}/boost_${BOOST_VERSION//\./_}.tar.gz"
BOOST_DIR="${DEPS_DIR}/boost-${BOOST_VERSION}"
BOOST_ROOT="/usr"
export BOOST_DIR
export BOOST_ROOT

# we cache this boost directory so we don't build it every time
if [[ -z "$(ls -A "${BOOST_DIR}")" ]]; then
  mkdir -p "${BOOST_DIR}"
  { wget --quiet -O - "${BOOST_URL}" | tar --strip-components=1 -xz -C "${BOOST_DIR}"; } || exit 1
fi

cd "$BOOST_DIR" || exit 1
{ sudo ./bootstrap.sh --with-libraries="$BOOST_LIBS" --prefix="$BOOST_ROOT" && sudo ./b2 install; } || exit 1
