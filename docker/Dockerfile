FROM ubuntu:24.04

LABEL org.opencontainers.image.created="2023-03-11"
LABEL org.opencontainers.image.url="https://github.com/odelaneau/shapeit5"
LABEL org.opencontainers.image.version="5.1.0"
LABEL org.opencontainers.image.licences="MIT"
LABEL org.opencontainers.image.title="shapeit5"
LABEL org.opencontainers.image.authors="olivier.delaneau@gmail.com"

# set the environment variables

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ARG DEBIAN_FRONTEND=noninteractive

RUN chsh -s /bin/bash
# run update and install necessary tools from package manager
WORKDIR /tmp
RUN apt -y update && apt -y -y upgrade && apt -y install -y \
	apt-utils bash time gnupg gcc g++ make autoconf \
	git zlib1g-dev bzip2 libbz2-dev liblzma-dev \
	build-essential git wget \
        libcurl4-gnutls-dev libgcrypt20-dev \
        libboost-all-dev libboost-iostreams-dev libhts-dev libhts3 && \
	apt -y clean && apt -y purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /user/src
COPY . shapeit5

WORKDIR /user/src/shapeit5

#Compile phase 1
RUN cd phase_common && make -j system && cp bin/phase_common ../static_bins/. && cd ..

#Compile phase 2
RUN cd phase_rare && make -j system && cp bin/phase_rare ../static_bins/. && cd ..

#Compile switch
RUN cd switch && make -j system && cp bin/switch ../static_bins/. && cd ..

#Compile ligate
RUN cd ligate && make -j system && cp bin/ligate ../static_bins/. && cd ..

#Compile simulate
RUN cd simulate && make -j system && cp bin/simulate ../static_bins/. && cd ..

#Compile xcftools
RUN cd xcftools && make -j system && cp bin/xcftools ../static_bins/.  && cd ..

RUN chmod +x /bin/*

