# PANDORA
# Pan-genome inference and genotyping with long noisy or short accurate reads

FROM ubuntu:bionic

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt update \
    && apt install -y software-properties-common \
    && apt-add-repository universe \
    && apt update \
    && apt install --no-install-recommends -y build-essential git cmake wget zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

#============================================
# INSTALL BOOST
#============================================
ENV BM 1
ENV Bm 62
ENV BP 0
ENV BOOST_VERSION "${BM}.${Bm}.${BP}"
ENV BOOST_V_USCORE "${BM}_${Bm}_${BP}"
ENV BOOST_URL "http://sourceforge.net/projects/boost/files/boost/${BOOST_VERSION}/boost_${BOOST_V_USCORE}.tar.gz"
ENV BOOST_LIBS "system,filesystem,iostreams,log,thread,date_time"
ENV BOOST_DIR "/boost_$BOOST_VERSION"
RUN mkdir -p "$BOOST_DIR" \
    && { wget --quiet -O - "${BOOST_URL}" | tar --strip-components=1 -xz -C "${BOOST_DIR}"; } \
    && apt-get remove -y wget
WORKDIR "$BOOST_DIR"
RUN ./bootstrap.sh --prefix=/usr/ --with-libraries="$BOOST_LIBS" \
    && ./b2 install \
    && cd .. \
    && rm -rf "$BOOST_DIR"
#============================================
# INSTALL PANDORA
#============================================
# can override the build type with docker's --build-arg command
# https://docs.docker.com/engine/reference/builder/#arg
ARG PANDORA_BUILD_TYPE="Release"
ENV PANDORA_DIR "/pandora/"

COPY . $PANDORA_DIR
WORKDIR ${PANDORA_DIR}/build
RUN cmake -DCMAKE_BUILD_TYPE="$PANDORA_BUILD_TYPE" .. \
    && make -j4 \
#    && ctest -V \
    && apt-get remove -y cmake git \
    && mv pandora /bin/pandora \
    && cd / \
    && rm -rf $PANDORA_DIR
