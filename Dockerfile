# PANDORA
# Pan-genome inference and genotyping with long noisy or short accurate reads

FROM ubuntu:20.04

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt update \
    && apt install -y software-properties-common \
    && apt-add-repository universe \
    && apt update \
    && apt install --no-install-recommends -y build-essential git cmake wget gdb \
                   minimap2=2.17+dfsg-2 mafft=7.453-1 racon=1.4.10-1build1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

#============================================
# INSTALL PANDORA
#============================================
# can override the build type with docker's --build-arg command
# https://docs.docker.com/engine/reference/builder/#arg
# ARG PANDORA_BUILD_TYPE="Release"
ARG PANDORA_BUILD_TYPE="Debug"
ENV PANDORA_DIR "/pandora/"

COPY . $PANDORA_DIR
WORKDIR ${PANDORA_DIR}/build
RUN cmake -DCMAKE_BUILD_TYPE="$PANDORA_BUILD_TYPE" -DHUNTER_JOBS_NUMBER=4 .. \
    && make -j4 \
    # && ctest -V \  # TODO: add tests back
    && apt-get remove -y cmake git \
    && mv pandora /bin/pandora \
    && cd / \
    && rm -rf $PANDORA_DIR \
    && rm -rf /root/.hunter

