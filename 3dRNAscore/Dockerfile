FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make

COPY 3dRNAscore /opt/3dRNAscore

RUN cd /opt/3dRNAscore \
 && make -j $(nproc)

ENV RNAscore=/opt/3dRNAscore \
    PATH="${PATH}:/opt/3dRNAscore/bin"
