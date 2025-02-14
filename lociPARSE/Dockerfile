# Use nvidia/cuda image
FROM nvidia/cuda:12.3.1-devel-ubuntu22.04
ENV TZ=Europe/Brussels
ENV DEBIAN_FRONTEND=noninteractive

# Set bash as current shell
RUN chsh -s /bin/bash
SHELL ["/bin/bash", "-c"]

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Anaconda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh -O ~/anaconda.sh \
    && /bin/bash ~/anaconda.sh -b -p /opt/conda \
    && rm ~/anaconda.sh \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy

# Set path to conda
ENV PATH /opt/conda/bin:$PATH

# Copy project files
COPY lociPARSE-main /opt/lociPARSE-main

# Create and activate conda environment
WORKDIR /opt/lociPARSE-main
RUN conda env create -f lociPARSE_environment.yml

# Set default command to activate conda environment
SHELL ["conda", "run", "-n", "lociPARSE", "/bin/bash", "-c"]
