
FROM ubuntu:22.04

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y

#/path/to/the/rasp-fd/folder
ENV RASP /home/rasp-fd-1.0/
RUN export RASP
COPY . /home/
RUN cd /home/rasp-fd-1.0/ && make
#RUN make
