FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make \
 && rm -rf /var/lib/apt/lists/*

COPY cgRNASP-CN opt/cgRNASP-CN
COPY bin/cgRNASP-CN /usr/bin/cgRNASP-CN

RUN gcc /opt/cgRNASP-CN/cgRNASP-CN.c -lm -o /opt/cgRNASP-CN/cgRNASP-CN
