## João Raimundo @BioData.pt

# Create a new environment to install R v4.0.0 (optional)_____________________//
 
conda create --name R_4.0.0  # create a new conda env for testing R v4.0.0
 
conda activate R_4.0.0
 
# install R v4.0.0 with conda________________________________________________//
 
 conda install -c conda-forge r-base=4.0.0 -y 
 
R --version # check R version (4.0.0)
 
R  # enter in R shell

# install DADA2_____________________________________________________________//
 
BiocManager::install(version = '3.11')
 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")

packageVersion("dada2") # # check DADA2 version (1.16.0)

# install dependencies

install.packages("ggplot2")
