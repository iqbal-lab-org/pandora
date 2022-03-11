# Based on https://github.com/radupopescu/musl-builder/blob/master/Dockerfile
# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora_static_binary_toolchain_musl:0.0.2

# This container has the C++ toolchain to build a static binary for pandora
# to build: docker build . -t leandroishilima/pandora_static_binary_toolchain_musl:0.0.2
FROM alpine:3.15

RUN apk update && apk add --no-cache binutils cmake make libgcc musl-dev gcc g++ bash linux-headers libexecinfo-dev git
