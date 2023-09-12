
FROM ubuntu:22.04

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y



COPY . /home/
ENV RNAscore /home/3dRNAscore/
RUN export RNAscore
RUN cd /home/3dRNAscore/ && make
#RUN make
