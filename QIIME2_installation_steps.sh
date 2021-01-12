#!/bin/bash
# João Raimundo @BioData.pt

###_Steps_to_QIIME2_________________________________________________________________________________________________//

##_Installing_with_Miniconda:___________________________________________________________________________________/

# Miniconda Installation:
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download miniconda
$ chmod +x Miniconda3-latest-Linux-x86_64.sh # define as executable
$ ./Miniconda3-latest-Linux-x86_64.sh -y # installation

# QIIME2 Installation with Conda
$ wget https://data.qiime2.org/distro/core/qiime2-2020.11-py36-linux-conda.yml
conda env create -n qiime2 --file qiime2-2020.11-py36-linux-conda.yml
# OPTIONAL CLEANUP
$ rm qiime2-2020.11-py36-linux-conda.yml

