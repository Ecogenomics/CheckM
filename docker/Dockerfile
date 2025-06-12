FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt install -y git python3-pip hmmer prodigal pplacer wget \
    && pip3 install numpy matplotlib pysam checkm-genome

RUN wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && tar -zxvf checkm_data_2015_01_16.tar.gz 
RUN rm checkm_data_2015_01_16.tar.gz

ENV CHECKM_DATA_PATH=/
ENV PATH="$HOME/.local/bin:$PATH"
CMD checkm	
 
