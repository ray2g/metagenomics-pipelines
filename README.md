# Metagenomic Pipelines

<br>

## DADA2

### Inputs:
* fastq samples (fastq)
* Silva_train_set.fa.gz (assigning taxonomy)
* Silva_species.fa.gz  (assigning taxonomy)

### Outputs:
* number_sequences_before_and_after_trimming.tsv (filter/trimming)
* asv_table.tsv (asv table with chimeras)
* summary_sequences_absolute.tsv (Track Reads through the pipeline -absolute)
* summary_sequences_percentage.tsv (Track Reads through the pipeline -percentage)
* asvTblNoChim.tsv (asv table without chimeras)
* taxTbl.tsv (taxonomy table)
* asvTaxTbl.tsv (asv + taxonomy table merged)
* asvTaxTbl.biom (asv + taxonomy table merged in BIOM format)





