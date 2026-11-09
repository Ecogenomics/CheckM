FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

# gPlace (PyPI: gtdb-gplace) replaced pplacer in CheckM 1.3.0 and requires Python >= 3.12
RUN apt-get update \
    && apt install -y git python3-pip hmmer prodigal wget \
    && pip3 install --break-system-packages numpy matplotlib pysam checkm-genome "gtdb-gplace>=0.1.0"

RUN wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && tar -zxvf checkm_data_2015_01_16.tar.gz 
RUN rm checkm_data_2015_01_16.tar.gz

ENV CHECKM_DATA_PATH=/
ENV PATH="$HOME/.local/bin:$PATH"
CMD checkm	
 
