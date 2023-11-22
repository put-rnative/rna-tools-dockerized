FROM ubuntu:22.04

RUN apt-get update \
 && apt-get install -y \
        build-essential \
        make \
 && rm -rf /var/lib/apt/lists/*

COPY cgRNASP /opt/cgRNASP
COPY bin/cgRNASP /usr/bin/cgRNASP
COPY bin/cgRNASP-C /usr/bin/cgRNASP-C
COPY bin/cgRNASP-PC /usr/bin/cgRNASP-PC

RUN gcc /opt/cgRNASP/cgRNASP/cgRNASP.c -lm -o /opt/cgRNASP/cgRNASP/cgRNASP \
 && gcc /opt/cgRNASP/cgRNASP-C/cgRNASP-C.c -lm -o /opt/cgRNASP/cgRNASP-C/cgRNASP-C \
 && gcc /opt/cgRNASP/cgRNASP-PC/cgRNASP-PC.c -lm -o /opt/cgRNASP/cgRNASP-PC/cgRNASP-PC
