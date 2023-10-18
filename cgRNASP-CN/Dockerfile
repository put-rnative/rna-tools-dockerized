FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make

COPY cgRNASP-CN opt/cgRNASP-CN

RUN gcc /opt/cgRNASP-CN/cgRNASP-CN.c -lm -o /usr/bin/cgRNASP-CN
