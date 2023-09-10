
FROM ubuntu:22.04

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y



COPY . /home/
RUN cd /home/cgRNASP-CN/ && gcc cgRNASP-CN.c -lm -o cgRNASP-CN
#RUN make
