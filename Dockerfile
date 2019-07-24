# PANDORA
# Pan-genome inference and genotyping with long noisy or short accurate reads

ARG UBUNTU_VERSION="18.04"
FROM ubuntu:"$UBUNTU_VERSION"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

WORKDIR /usr/local

RUN apt update && apt install -y software-properties-common
RUN apt-add-repository universe && apt update
RUN apt install -y \
    build-essential \
    cmake \
    git \
    man \
    seqtk \
    time \
    wget


#============================================
# INSTALL ZLIB
#============================================
ENV ZLIB_VERSION 1.2.11
ENV ZLIB_URL="http://www.zlib.net/zlib-${ZLIB_VERSION}.tar.gz"

RUN wget "$ZLIB_URL" -O - | tar xzf -
WORKDIR zlib-"$ZLIB_VERSION"
RUN ./configure --prefix=/usr/ && make && make install

WORKDIR /usr/local

#============================================
# INSTALL BOOST
#============================================
ENV BOOST_MAJOR 1
ENV BOOST_MINOR 62
ENV BOOST_PATCH 0
ENV BOOST_URL="https://sourceforge.net/projects/boost/files/boost/${BOOST_MAJOR}.${BOOST_MINOR}.${BOOST_PATCH}/boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}.tar.gz"
ENV BOOST_LIBS="system,filesystem,iostreams,log,thread,date_time"

RUN wget "$BOOST_URL" -O - | tar xzf -
WORKDIR "boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}"
RUN ./bootstrap.sh --prefix=/usr/ --with-libraries="$BOOST_LIBS" && ./b2 install

WORKDIR /usr/local

#============================================
# INSTALL PANDORA
#============================================
ENV PANDORA_BRANCH dev
ENV PANDORA_GIT="https://github.com/rmcolq/pandora.git"

RUN git clone -b "$PANDORA_BRANCH" --single-branch  --recursive "$PANDORA_GIT"
WORKDIR pandora
RUN mkdir -p build
WORKDIR build
RUN cmake -DCMAKE_BUILD_TYPE=Release .. && make
RUN ln -s $(realpath pandora) /usr/local/bin/pandora
# uncomment the below once tests pass for release mode
# RUN ctest -V
