
FROM ubuntu:22.04

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y

RUN apt-get update \
  && apt-get -y install build-essential \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget https://github.com/Kitware/CMake/releases/download/v3.24.1/cmake-3.24.1-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod u+x /tmp/cmake-install.sh \
      && mkdir /opt/cmake-3.24.1 \
      && /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake-3.24.1 \
      && rm /tmp/cmake-install.sh \
      && ln -s /opt/cmake-3.24.1/bin/* /usr/local/bin


COPY . /home/
ENV BRiQ_DATAPATH=/home/RNA-BRiQ-data
ENV PATH=$PATH:/home/RNA-BRiQ//build/bin

#readme advises:
#wget http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz
#export BRiQ_DATAPATH=FILEPATH/BRiQ_data     ## Change "FILEPATH" to the real path to contain the preprocessed data
#tar -xvz RNA-BRiQ-data.tar.gz --one-top-level=$BRiQ_DATAPATH

# 
RUN cd /home/ && wget https://apisz.sparks-lab.org:8443/downloads/Resource/RNA/2_RNA_structure_prediction/RNA-BRiQ-data.tar.gz
RUN cd /home/RNA-BRiQ/ && mkdir build && cd build && cmake ../ && make \
	&& tar -xvf /home/RNA-BRiQ-data.tar.gz /home/RNA-BRiQ-data \   
	&& export PATH=$PATH:/home/RNA-BRiQ//build/bin 
#RUN make
