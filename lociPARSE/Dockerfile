# Use nvidia/cuda image
#https://developer.download.nvidia.com/licenses/NVIDIA_Deep_Learning_Container_License.pdf <- licencja
FROM nvidia/cuda:11.3.1-devel-ubuntu20.04
ENV TZ=Europe/Brussels

# set bash as current shell
RUN chsh -s /bin/bash
SHELL ["/bin/bash", "-c"]

# install anaconda
RUN apt-get update
RUN apt-get install -y wget bzip2 ca-certificates && apt-get clean
#libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion && \
#        apt-get clean
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

# set path to conda
ENV PATH /opt/conda/bin:$PATH


# setup conda virtual environment



#FROM continuumio/miniconda3:latest
#ENV DEBIAN_FRONTEND noninteractive

#RUN conda upgrade conda -y

COPY lociPARSE-main /opt/lociPARSE-main


#RUN cd /opt/lociPARSE-main && conda env create -f lociPARSE_environment.yml
#RUN conda activate lociPARSE