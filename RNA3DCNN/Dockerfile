



FROM continuumio/anaconda3:latest
RUN apt-get update && apt-get install -y build-essential libgl1-mesa-glx && \
    rm -rf /var/lib/apt/lists/*
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility

RUN apt-get update && \
      apt-get -y install sudo
	  
RUN apt update && apt install -y make
RUN apt-get update && apt-get install build-essential -y

# Install some dependencies
RUN apt-get update && apt-get install -y \
		bc \
		build-essential \
		cmake \
		curl \
		g++ \
		&& \
	apt-get clean && \
	apt-get autoremove && \
	rm -rf /var/lib/apt/lists/* 

# Install pip
RUN curl -O https://bootstrap.pypa.io/get-pip.py && \
	python get-pip.py && \
	rm get-pip.py

# Add SNI support to Python
RUN pip --no-cache-dir install \
		pyopenssl \
		ndg-httpsclient \
		pyasn1


# Install other useful Python packages using pip
RUN pip --no-cache-dir install --upgrade ipython && \
	pip --no-cache-dir install \
		jupyter \
		biopython\
		&& \
	python -m ipykernel.kernelspec

RUN pip --no-cache-dir install \
		tensorflow
	
RUN pip --no-cache-dir install \
		numpy
		
RUN pip --no-cache-dir install \
		keras
RUN conda install biopython 

# Install Keras
#RUN pip --no-cache-dir install git+git://github.com/fchollet/keras.git
#@${KERAS_VERSION}
COPY . /home/

