
FROM ubuntu:22.04

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y



COPY . /home/
RUN cd /home/cgRNASP/cgRNASP-C && gcc cgRNASP-C.c -lm -o cgRNASP-C
RUN cd /home/cgRNASP/cgRNASP-PC && gcc cgRNASP-PC.c -lm -o cgRNASP-PC
RUN cd /home/cgRNASP/cgRNASP &&  gcc cgRNASP.c -lm -o cgRNASP
#RUN make
