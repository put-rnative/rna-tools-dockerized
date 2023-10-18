FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make \
 && rm -rf /var/lib/apt/lists/*

COPY rsRNASP /opt/rsRNASP

RUN gcc /opt/rsRNASP/rsRNASP.c -lm -o /usr/bin/rsRNASP
