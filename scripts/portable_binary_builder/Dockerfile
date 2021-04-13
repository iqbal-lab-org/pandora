# based on https://github.com/radupopescu/musl-builder/blob/master/Dockerfile
# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora_static_binary_toolchain:0.0.1

# This container has the C++ toolchain to build a static binary for pandora
# to build: sudo docker build . -t leandroishilima/pandora_static_binary_toolchain:0.0.1
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y binutils-dev cmake make musl-dev gcc g++

