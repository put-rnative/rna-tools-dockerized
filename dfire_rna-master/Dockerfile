
FROM ubuntu:22.04

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y


COPY . /home/
RUN cd /home/dfire/ && make #&& install.sh
#RUN /home/dfire/install.sh
ENV DFIRE_RNA_HOME='/home/dfire/'
ENV foo /bar
# or ENV foo=/bar
#RUN make
