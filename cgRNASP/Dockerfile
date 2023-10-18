FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make \
 && rm -rf /var/lib/apt/lists/*

COPY cgRNASP /opt/cgRNASP

RUN gcc /opt/cgRNASP/cgRNASP/cgRNASP.c -lm -o /usr/bin/cgRNASP \
 && gcc /opt/cgRNASP/cgRNASP-C/cgRNASP-C.c -lm -o /usr/bin/cgRNASP-C \
 && gcc /opt/cgRNASP/cgRNASP-PC/cgRNASP-PC.c -lm -o /usr/bin/cgRNASP-PC
