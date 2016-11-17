FROM bioconductor/release_base
MAINTAINER "Christian Diener <mail[at]cdiener.com>"

RUN apt-get update && apt-get -y -t unstable install python3 python3-pip libssl-dev wget && apt-get clean
RUN pip3 install python-libsbml numpy scipy lxml cobra pandas && rm -rf /tmp/*

RUN mkdir /data && chown rstudio:rstudio /data

USER rstudio

RUN git clone https://github.com/cdiener/proliferation /data/proliferation
WORKDIR /data/proliferation
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); \
    biocLite('BiocInstaller'); setRepositories(ind=1:2); \
    install.packages(c('devtools', 'rmarkdown')); \
    devtools::install_deps('prtools', dependencies=TRUE); \
    devtools::install('prtools')"
RUN wget https://zenodo.org/record/61980/files/regprob.csv && \
    wget https://zenodo.org/record/61982/files/tcga.rds

USER root
