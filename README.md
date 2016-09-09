[![wercker status](https://app.wercker.com/status/4c8247e9636b875cb647a4173200b674/s "wercker status")](https://app.wercker.com/project/bykey/4c8247e9636b875cb647a4173200b674)
![paper status](https://img.shields.io/badge/paper-in_preparation-yellow.svg)
[![protocol status](https://img.shields.io/badge/protocol-online-green.svg)](https://cdiener.github.io/proliferation)

# Pan-cancer analysis of proliferation rates and metabolic liabilities

Project to study proliferation and metabolic liabilities across several
thousand cancer samples. The automatically generated output can be found
at https://cdiener.github.io/proliferation

## What this repository contains

This repository contains all additional files to fully reproduce all analyses.
It uses a [R markdown](http://rmarkdown.rstudio.com/) file to automatically
create a detailed protocol of the analysis which is available at
https://cdiener.github.io/proliferation.

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

For Debian Jessie we recommend installation of Python via [Anaconda](https://www.continuum.io/downloads)
in order to get a version of pip that is greater than 8.1.

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
[Amazon AWS](https://aws.amazon.com/) you can use the docker image.

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
5. Access the machine at http://<your-ip>:8000 where <your-ip> is the IP of
   your VM or "localhost" when running docker on your own machine. You will
   be prompted with for login information where the user and password are
   "rstudio".

This will present you with an R studio interface where all additional dependencies
and intermediate data are available.

For a start select the folder "proliferation" in the file panel on the lower
right and click on "protocol.rmd".
