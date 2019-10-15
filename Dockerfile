FROM continuumio/miniconda3:4.6.14

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython python3-pip wget

RUN pip3 install biopython ete3 pycosat PyYAML requests

# Install bowtie2 and krakenuniq
RUN conda install -c bioconda bowtie2 && \
    conda install -c bioconda krakenuniq && \
    conda clean -afy
