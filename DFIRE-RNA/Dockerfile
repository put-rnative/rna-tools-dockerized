FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make

COPY dfire /opt/dfire

RUN cd /opt/dfire \
 && make -j $(nproc)

ENV DFIRE_RNA_HOME=/opt/dfire \
    PATH="${PATH}:/opt/dfire/bin"
