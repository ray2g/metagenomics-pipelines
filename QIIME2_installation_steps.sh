#!/bin/bash
# João Raimundo @BioData.pt

###_Steps_to_QIIME2_________________________________________________________________________________________________//

##_Installing_with_Miniconda:___________________________________________________________________________________/

# Miniconda Installation:
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download miniconda
$ chmod +x Miniconda3-latest-Linux-x86_64.sh # define as executable
$ ./Miniconda3-latest-Linux-x86_64.sh -y # installation

# Create a new conda environment
$ sudo apt-get update && apt-get upgrade -y
$ conda create --name qiime2-2020.11 -y # in this example we choose "qiime_2020.11" as the name.
$ conda activate qiime2-2020.11 # activate the environment

# QIIME2 Installation
$ conda install -c qiime2/label/r2020.11 qiime2 -y # installing qiime2 release 2020.11 in the new environment

