FROM ubuntu:20.04
MAINTAINER Karthik G <gkarthik@scripps.edu>

RUN apt-get update
RUN apt-get install -y build-essential autoconf zlib1g-dev python3 wget libbz2-dev liblzma-dev libncurses-dev git bedtools python3-pip vim nano
# HTSlib
RUN cd root/ &&\
    wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 &&\
    tar xvf htslib-1.10.2.tar.bz2 &&\
    cd htslib-1.10.2/ &&\
    ./configure &&\
    make &&\
    make install &&\
    cd ../ &&\
    rm htslib-1.10.2.tar.bz2
ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
# SAMtools
RUN cd root &&\
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 &&\
    tar xvf samtools-1.9.tar.bz2 &&\
    cd samtools-1.9/ &&\
    ./configure &&\
    make &&\
    make install &&\
    cd ../ &&\
    rm samtools-1.9.tar.bz2
# iVar
RUN cd root/ &&\
    git clone https://github.com/andersen-lab/ivar.git &&\
    cd ivar/ &&\
    ./autogen.sh &&\
    ./configure &&\
    make &&\
    make install
# bwa
RUN cd root/ &&\
    wget https://github.com/lh3/bwa/archive/v0.7.18.tar.gz &&\
    tar xvf v0.7.18.tar.gz &&\
    cd bwa-0.7.18/ &&\
    make &&\
    cd ../ &&\
    rm v0.7.18.tar.gz
ENV PATH /root/bwa-0.7.18:$PATH
# Snakemake
RUN pip3 install pandas snakemake
