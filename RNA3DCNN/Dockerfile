FROM continuumio/miniconda3:latest

RUN conda upgrade conda \
 && conda install -c conda-forge \
        biopython \
        keras \
        tensorflow

COPY RNA3DCNN /opt/RNA3DCNN
COPY rna3dcnn /usr/bin/rna3dcnn
