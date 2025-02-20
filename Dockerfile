FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get install -y \
      build-essential \
      cmake \
      make \
      parallel \
      python3-biopython \
      python3-pandas \
      python3-pip \
      python3-setuptools \
      python3-virtualenv \
 && rm -rf /var/lib/apt/lists/*

# -O2: turn on optimization level 2 is commonly used for production—and has been extensively tested
# -march=native: target native architecture
# -flto: enable Link Time Optimization
# -fstack-protector-strong: inserts guards on functions that use potentially dangerous constructs, reducing the risk of stack smashing exploits
# -D_FORTIFY_SOURCE=2: when used with optimization, this macro helps detect some buffer overflows and format–string bugs in certain glibc functions
# -fPIE: tells the compiler to produce position–independent code, a prerequisite for enabling address space layout randomization (ASLR)
# -fstack-clash-protection: helps mitigate stack clash vulnerabilities
# -Wl,-z,relro: “Read–only relocations” helps protect memory regions by marking them read–only after initialization
# -Wl,-z,now: forces immediate binding of symbols (i.e. no lazy binding) - reducing the risk of certain types of runtime attacks
# -pie: produces an executable that is position–independent, allowing the operating system to randomize its location in memory (ASLR)
# -Wl,--as-needed: causes the linker to only add libraries that are actually used, reducing the code footprint
ENV CFLAGS="-O2 -march=native -flto -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIE -fstack-clash-protection" \
    CXXFLAGS="-O2 -march=native -flto -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIE -fstack-clash-protection" \
    LDFLAGS="-flto -Wl,-z,relro -Wl,-z,now -pie -Wl,--as-needed" \
    RNAscore=/opt/3dRNAscore \
    BRiQ_DATAPATH=/opt/BRiQ_data \
    DFIRE_RNA_HOME=/opt/dfire \
    RASP=/opt/rasp-fd-1.0


COPY --link \
  3dRNAscore/ \
  BRiQ/ \
  DFIRE-RNA/ \
  RASP/ \
  RNA3DCNN/ \
  cgRNASP-CN/ \
  cgRNASP/ \
  lociPARSE/ \
  rsRNASP/ \
  /opt/
ADD --link \
  BRiQ/RNA-BRiQ-data-atomicEne.tar.xz \
  BRiQ/RNA-BRiQ-data-baseAtomClash.tar.xz \
  BRiQ/RNA-BRiQ-data-split6D.tar.xz \
  BRiQ/RNA-BRiQ-data-other.tar.xz \
  /opt/


# Compiled
RUN echo "=== 3dRNAscore ===" \
 && cd /opt/3dRNAscore \
 && make -j $(nproc) \
 && echo "=== RNA-BRiQ ===" \
 && mkdir /opt/RNA-BRiQ/build \
 && cd /opt/RNA-BRiQ/build \
 && cmake .. \
 && make -j $(nproc) \
 && echo "=== cgRNASP, cgRNASP-CN, rsRNASP ===" \
 && gcc /opt/cgRNASP/cgRNASP/cgRNASP.c -lm -o /opt/cgRNASP/cgRNASP/cgRNASP \
 && gcc /opt/cgRNASP/cgRNASP-C/cgRNASP-C.c -lm -o /opt/cgRNASP/cgRNASP-C/cgRNASP-C \
 && gcc /opt/cgRNASP/cgRNASP-PC/cgRNASP-PC.c -lm -o /opt/cgRNASP/cgRNASP-PC/cgRNASP-PC \
 && gcc /opt/cgRNASP-CN/cgRNASP-CN.c -lm -o /opt/cgRNASP-CN/cgRNASP-CN \
 && gcc /opt/rsRNASP/rsRNASP.c -lm -o /opt/rsRNASP/rsRNASP \
 && echo "=== DFIRE ===" \
 && cd /opt/dfire \
 && make -j $(nproc) \
 && echo "=== RASP ===" \
 && cd /opt/rasp-fd-1.0 \
 && make -j $(nproc)


# Python-based
RUN echo "=== RNA3DCNN ===" \
 && python3 -m virtualenv /opt/RNA3DCNN/venv \
 && /opt/RNA3DCNN/venv/bin/pip install 'tensorflow<=2.15' \
 && echo "=== lociPARSE ===" \
 && cd /opt/lociPARSE-8c7acbe4e7c486122a4c261b1ea68fff7247b796 \
 && pip install .


# Wrapper
COPY scoring-wrapper.py /usr/bin/
