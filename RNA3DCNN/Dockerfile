FROM continuumio/miniconda3:24.1.2-0

RUN conda upgrade conda \
 && conda install -c conda-forge \
        biopython \
        'keras<3' \
        tensorflow

COPY RNA3DCNN /opt/RNA3DCNN
