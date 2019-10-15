FROM ubuntu:18.04

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython python3-pip wget

RUN pip3 install biopython ete3 pycosat PyYAML requests

RUN cd /usr/local/ && \
    mkdir miniconda && \
    cd miniconda && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    /bin/bash -c "source ~/.bashrc"

# Install bowtie2 and krakenuniq
RUN conda install -c bioconda bowtie2 && \
    conda install -c bioconda krakenuniq && \
    conda clean -afy
