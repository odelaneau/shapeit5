FROM ubuntu:20.04

LABEL org.opencontainers.image.created="2022-02-02"
LABEL org.opencontainers.image.url="https://github.com/odelaneau/GLIMPSE"
LABEL org.opencontainers.image.version="0.0.1"
LABEL org.opencontainers.image.licences="MIT"
LABEL org.opencontainers.image.title="glimpse2"
LABEL org.opencontainers.image.authors="simone.rubinacci@unil.ch"

# set the environment variables
ENV HOME /root

ENV htslib_version 1.14
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ARG DEBIAN_FRONTEND=noninteractive

# run update and install necessary tools from package manager
WORKDIR /tmp
RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	apt-utils curl gnupg gcc g++ make autoconf && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	git zlib1g-dev bzip2 libbz2-dev liblzma-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# bcftools
RUN curl -SL https://github.com/samtools/bcftools/releases/download/${htslib_version}/bcftools-${htslib_version}.tar.bz2 \
    > bcftools-${htslib_version}.tar.bz2 && \
    tar -xjvf bcftools-${htslib_version}.tar.bz2 && \
    cd /tmp/bcftools-${htslib_version} && \
    ./configure && \
    make && \
    make install && \
    rm -rf bcftools-${htslib_version}.tar.bz2 && \
    rm -rf bcftools-${htslib_version}

# get GLIMPSE v2.0.0 [binary]
ADD resources/GLIMPSE_phase_static /bin/GLIMPSE_phase_v2.0.0
ADD resources/GLIMPSE_split_reference_static /bin/GLIMPSE_split_reference_v2.0.0
ADD resources/GLIMPSE_call_static /bin/GLIMPSE_call_v2.0.0
ADD resources/GLIMPSE_concordance_static /bin/GLIMPSE_concordance_v2.0.0

RUN chmod +x /bin/GLIMPSE*

WORKDIR $HOME
