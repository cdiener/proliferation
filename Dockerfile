FROM bioconductor/release_base
MAINTAINER "Christian Diener <mail[at]cdiener.com>"

RUN apt-get -y -t testing python3 python3-pip git wget && apt-get clean
RUN pip3 install python-libsbml numpy scipy lxml cobra && rm -rf /tmp/*

USER rstudio

RUN cd /home/rstudio
RUN git clone https://github.com/cdiener/proliferation && cd proliferation
RUN RUN Rscript -e "install.packages('devtools');
    devtools::install_deps('prtools', dependencies=TRUE);
    devtools::install('prtools')"
RUN wget https://zenodo.org/record/61980/files/regprob.csv && \
    wget https://zenodo.org/record/61982/files/tcga.rds
