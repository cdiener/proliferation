[![wercker status](https://app.wercker.com/status/4c8247e9636b875cb647a4173200b674/s "wercker status")](https://app.wercker.com/project/bykey/4c8247e9636b875cb647a4173200b674)
[![codecov](https://codecov.io/gh/cdiener/proliferation/branch/master/graph/badge.svg)](https://codecov.io/gh/cdiener/proliferation)
![paper status](https://img.shields.io/badge/paper-submitted-yellow.svg)
[![Docker Pulls](https://img.shields.io/docker/pulls/cdiener/proliferation.svg?maxAge=2592000)](https://hub.docker.com/r/cdiener/proliferation)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.154546.svg)](https://doi.org/10.5281/zenodo.154546)


# Pan-cancer analysis of proliferation rates and metabolic liabilities

Project to study proliferation and metabolic liabilities across several
thousand cancer samples. The automatically generated output can be found
at https://cdiener.github.io/proliferation

## What this repository contains

This repository contains all additional files to fully reproduce all analyses.
It uses a [R markdown](http://rmarkdown.rstudio.com/) file to automatically
create a detailed protocol of the analysis which is available at
https://cdiener.github.io/proliferation.

## Intermediate data

In case you would like to reproduce the analysis but do now want to wait
downloading the raw data we also provide some intermediate data files you
can (optinonally) use to accelerate the analysis:

1. NCI-60 gene expression data and proliferation rates
   http://dx.doi.org/10.5281/zenodo.61980
2. TCGA data as compressed RDS file
   http://dx.doi.org/10.5281/zenodo.61982
3. probe map as compressed RDS file
   *Already in the Github repository*

## Proliferation rate predictions

The predicted proliferation rates for more than 12000 samples from the TCGA.

1. Raw predictions for 12339 samples together with their respective barcode and cancer panel 
   https://github.com/cdiener/proliferation/results/pred_rates.csv
2. Raw predictions for 12111 samples that had clinical annotations merged with several clinical
   indicators such as survival data, histology types, TNM stage, follow ups and no. of new tumor
   events.
   https://github.com/cdiener/proliferation/results/combined.csv

## Rerunning the analysis interactively

### Locally

**Note that this will require a machine with 16+ GB of RAM.**

If you have docker installed you can use the steps provided in the section below
"With a cloud provider".

For a local installation you will need R (http://r-project.org), git
(http://git-scm.org) and Python 3 installed (http://python.org).

On Ubuntu or Debian Stretch these can be installed via

```bash
sudo apt-get install r-base r-base-dev python3 python3-pip git
```

Clone the repository and enter the folder:
```bash
git clone https://github.com/cdiener/proliferation && cd proliferation
```

For Debian Jessie we recommend updating pip on a per-user setting
in order to get a version of pip that is greater than 8.1.

```bash
sudo apt-get install r-base r-base-dev python3 python3-pip git
pip3 install -U pip
```

The python dependencies can be installed with

```bash
pip3 install lxml python-libsbml numpy scipy cobra
```

R dependencies can be installed from within R (type "R" in Terminal) with

```R
install.packages("devtools")
devtools::install_github("cdiener/proliferation/prtools")
```

You can now run the steps presented in the protocol.

### With a cloud provider

Using a cloud provider such as [Google Cloud](https://cloud.google.com/) or
[Amazon AWS](https://aws.amazon.com/) you can rerun the analysis easily with
the provided docker image.

1. Create a new virtual machine with more than 16 GB of RAM using the CoreOS
   stable image.
2. Login to the machine using SSH as described by your cloud provider.
3. Get the docker image with

   ```bash
   docker pull cdiener/proliferation
   ```
4. Run the docker image

   ```bash
   docker run -d -p 8000:8787 cdiener/proliferation
   ```

   (You can now disconnect the SSH connection, but leave the VM running.)
5. Access the machine at http://your-ip:8000 where "your-ip" is the IP of
   your VM or "localhost" when running docker on your own machine. You will
   be prompted with for login information where the user and password are
   "rstudio".

This will present you with an R studio interface where all additional dependencies
and intermediate data are available.

For a start click the "..." symbol in the file panel on the lower
right and enter "/data/proliferation". Now click on "protocol.rmd" to see or
run the protocol (by clicking "knit HTML") or use any of the *.R
files in the same directory.
