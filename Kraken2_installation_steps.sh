#!/bin/bash
###_Steps_to_install_Kraken_2_________________________________________________________________________//

##_Installing_with_Miniconda:______________________________________________________________________/

# Miniconda Installation:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download miniconda
chmod +x Miniconda3-latest-Linux-x86_64.sh # set as executable
./Miniconda3-latest-Linux-x86_64.sh -y # installation

# Create a new conda environment
sudo apt-get update && apt-get upgrade -y
conda create --name kraken2 -y # in this example we choose "kraken2" as the name
conda activate kraken2 # activate the environment

# Kraken 2 Installation
conda install -c bioconda kraken2 -y # installing kraken2 in the new environment

