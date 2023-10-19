FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        cmake \
        make \
 && rm -rf /var/lib/apt/lists/*

COPY RNA-BRiQ /opt/RNA-BRiQ
ADD RNA-BRiQ-data-atomicEne.tar.xz /opt
ADD RNA-BRiQ-data-baseAtomClash.tar.xz /opt
ADD RNA-BRiQ-data-split6D.tar.xz /opt
ADD RNA-BRiQ-data-other.tar.xz /opt

RUN mkdir /opt/RNA-BRiQ/build \
 && cd /opt/RNA-BRiQ/build \
 && cmake .. \
 && make -j $(nproc)

ENV BRiQ_DATAPATH=/opt/BRiQ_data \
    PATH="${PATH}:/opt/RNA-BRiQ/build/bin"
