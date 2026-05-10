FROM continuumio/miniconda3

RUN conda install -y \
    -c bioconda \
    -c conda-forge \
    sra-tools \
    fastqc \
    fastp \
    spades \
    bwa \
    samtools \
    bcftools

WORKDIR /data