#!bin/bash
## João Raimundo @BioData.pt

# directory setup for DADA2 pipeline
mkdir <analysisName> <analysisName>/tables <analysisName>/plots <analysisName>/data \
<analysisName>/taxonomy <analysisName>/taxonomy/database <analysisName>/taxonomy/objects

# download SILVA v138 database: train set & species
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz?download=1 -O <analysisName>/taxonomy/database/silva_nr99_v138_train_set.fa.gz
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz?download=1 -O <analysisName>/taxonomy/database/silva_species_assignment_v138.fa.gz


