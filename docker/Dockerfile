
ARG base_image
FROM $base_image

ARG HTSLIB_VER=1.14


RUN apt-get update

RUN apt-get install wget make curl gcc g++ cmake gfortran libreadline6-dev libboost-all-dev libopenblas-dev default-jre lsb-core libcairo2-dev libpcre2-dev libpcre3-dev libssl-dev libcurl4-openssl-dev libz-dev libbz2-dev liblzma-dev zlib1g-dev --yes


RUN wget https://cloud.r-project.org/src/base/R-4/R-4.1.2.tar.gz && tar xvzf R-4.1.2.tar.gz \
  && cd R-4.1.2 && ./configure --with-x=no --enable-R-shlib=yes --with-blas="-lopenblas" && \
    make && \
    mkdir -p /usr/local/lib/R/lib && \
    make install && \
    cd .. && \
    rm -rf R-4.1.2*


RUN curl -LO https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 && \
    tar -xvjf htslib-${HTSLIB_VER}.tar.bz2 && cd htslib-${HTSLIB_VER}&& \
    ./configure && make && make install && cd .. && rm -rf htslib-${HTSLIB_VER}*


RUN apt-get install moreutils  zip unzip -y
RUN pip3 install --no-cache numpy scipy pandas cget
RUN pip3 install csvkit


COPY requirements.txt .
RUN apt-get install python3-pip -y && pip3 install -r requirements.txt


COPY scripts/qqplot.R /usr/local/bin
RUN chmod a+x /usr/local/bin/qqplot.R


COPY scripts/install.packages.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/install.packages.R && \
    install.packages.R

ADD ./scripts ../scripts
