FROM biocontainers/biocontainers:latest

RUN conda install -c bioconda bwa=0.7.17
RUN conda install -c bioconda samtools=1.5
RUN conda install -c bioconda picard=2.17.11
RUN conda install -c bioconda bamtools=2.4.0
RUN conda install -c bioconda ngs-disambiguate=2018.05.03
RUN conda install -c bioconda star=2.6.1a

WORKDIR /data
