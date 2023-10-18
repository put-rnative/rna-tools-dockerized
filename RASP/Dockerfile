FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make \
 && rm -rf /var/lib/apt/lists/*

COPY rasp-fd-1.0 /opt/rasp-fd-1.0

RUN cd /opt/rasp-fd-1.0 \
 && make

ENV RASP=/opt/rasp-fd-1.0 \
    PATH="${PATH}:/opt/rasp-fd-1.0/bin"
