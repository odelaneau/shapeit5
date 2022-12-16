FROM amazonlinux:2022 AS builder

#dependencies
RUN yum update -y && yum install -y gcc zlib-devel make tar bzip2 ncurses-devel bzip2-devel xz-devel gzip gcc-c++ libcurl-devel openssl-devel cmake boost-devel git boost-static findutils

RUN curl -O -L https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
RUN tar -xjvf bcftools-1.16.tar.bz2
RUN cd bcftools-1.16/htslib-1.16/ && ./configure --enable-s3 && make -j4 && make install && cd -
RUN cd bcftools-1.16/ && ./configure --enable-s3 && cd -
RUN make -C bcftools-1.16/ -j4 install

ADD . /shapeit5
RUN for dir in phase_common phase_rare switch ligate;do make -C /shapeit5/${dir} amazonlinux2022 -j4; done

FROM amazonlinux:2022

COPY --from=builder /usr/local/bin/bcftools /usr/bin/
COPY --from=builder /usr/local/bin/tabix /usr/bin/
COPY --from=builder /usr/local/bin/bgzip /usr/bin/
COPY --from=builder /shapeit5/*/bin/* /usr/bin/
