## Joao Raimundo @BioData.pt ----------------------///

#______Import_Packages________________________//

message('>>> Importing Packages... \n')

library("dada2") 
library("ggplot2")
library("DT")


#______Directories______________________________________________________________________________________________//

setwd("~/Desktop/dada2_test/") # set the working directory ## CHANGE ME

path <- "~/Desktop/dada2_test/data/" # set the data directory ## CHANGE ME  

plot_dir <- "~/Desktop/dada2_test/plots" # set the directory to save the plots ## CHANGE ME

fastqPath <- "~/Desktop/dada2_test/data/fastq" # set the directory to save processed fastq files ## CHANGE ME

tables_dir <- "~/Desktop/dada2_test/tables" # set the directory to save the tables ## CHANGE ME

taxonomy_dir <- "~/Desktop/dada2_test/taxonomy" # set the directory for taxonomy databases and to save R objects ## CHANGE ME

message("\n>> Reading data... \n") 
list.files(path)


#____Match_Forward_and_Reverse_fastq_files_____________________________________________//

#set the relative path to to each forward and reverse fastq files directory
fastqFwdPath <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fastqRevPath <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sampleNames <- sapply(strsplit(basename(fastqFwdPath), "_"), `[`, 1)


#____Inspect_read_quality_profiles________________________________________________________________________________//

message("\n>>> Inspecting Read Quality Profiles... \n")

fwd_qualityProfile <- plotQualityProfile(fl = fastqFwdPath) # forward 
rev_qualityProfile <- plotQualityProfile(fl = fastqRevPath) # reverse

message("\n>> Saving read quality profile plots on: ", paste0("\n"), plot_dir, paste0("/fwd_qualityProfile.png"), 
                                                     paste0("\n"), plot_dir, paste0("rev_qualityProfile.png"), paste0("\n"))

# saving quality profiles plots
ggsave(filename = "fwd_qualityProfile.png", plot = fwd_qualityProfile, path = plot_dir, width = 12, height = 8)
ggsave(filename = "rev_qualityProfile.png", plot = rev_qualityProfile , path = plot_dir, width = 12, height = 8)


#___Filter_and_trim_reads______________________________________________________________________________________________________________//

message(">>> Filter and Trimming the Reads... \n")
message("Forward Reads truncLen: 240")
message("Reverse Reads truncLen: 160 \n")

# relative file path for filtered Fwd and Rev reads that will be created below 
filtFastqFwdPath <- file.path(fastqPath, "filtered", paste0(sampleNames, "_fwd_filt.fastq.gz")) # Fwd
filtFastqRevPath <- file.path(fastqPath, "filtered", paste0(sampleNames, "_rev_filt.fastq.gz")) # Rev

# assign to each file path the sample name
names(filtFastqFwdPath) <- sampleNames # Fwd
names(filtFastqRevPath) <- sampleNames # Rev

# filter and trim fwd and rev fastq files writing the new filtered files in compressed
filterTrimReads <- filterAndTrim(fwd = fastqFwdPath, filt = filtFastqFwdPath, rev = fastqRevPath, filt.rev = filtFastqRevPath, 
                                 truncLen = c(240,160), maxEE = c(2,2), truncQ = 2, maxN = 0, rm.phix = TRUE, 
                                 compress = TRUE, verbose = FALSE, multithread = TRUE)

message(paste0(">> Filtered reads were saved on: \n"),fastqPath,paste0("/filtered/\n"))

message(">> Number of sequences kept after filtering/trimming in relation to the initial number of reads: \n")
filterTrimReads

message("\n")

message(">>> Inspecting Trimmed Read Quality Profiles... \n")

# forward quality plots after filtering and trimming
fwd_trimmed_qualityProfile <- plotQualityProfile(fl = filtFastqFwdPath)

# reverse quality plots after filtering and trimming
rev_trimmed_qualityProfile <- plotQualityProfile(fl = filtFastqRevPath)

message(">> Saving trimmed reads quality profile plots on: ", 
        paste0("\n"), plot_dir, paste0("/fwd_qualityProfile_trimmed.png"), 
        paste0("\n"), plot_dir, paste0("/rev_qualityProfile_trimmed.png"), paste0("\n"))

# saving trimmed reads quality profiles plots
ggsave(filename = "fwd_qualityProfile_trimmed.png", plot = fwd_trimmed_qualityProfile, path = plot_dir, width = 12, height = 8)
ggsave(filename = "rev_qualityProfile_trimmed.png", plot = rev_trimmed_qualityProfile , path = plot_dir, width = 12, height = 8)

message(">> Saving the trimmed summary table on: \n", 
        tables_dir, paste0("/number_sequences_before_and_after_trimming.csv\n"))

# save table filterTrimReads
write.csv(filterTrimReads,paste0(tables_dir,"/number_sequences_before_and_after_trimming.csv"))

#____Learning_the_Error_Rates_________________________________________________________________________________________//

message(">>> Learning the Error Rates... \n")

# model/learn the error rates for the filtered fastq files
message(">> Forward reads:")
errFwd <- learnErrors(fls = filtFastqFwdPath, multithread = TRUE) # fwd
message("\n>> Reversed reads:")
errRev <- learnErrors(fls = filtFastqRevPath, multithread = TRUE) # rev

# Plot errors 
fwd_plotErrors <- plotErrors(dq = errFwd, nominalQ = TRUE) # fwd
rev_plotErrors <- plotErrors(dq = errRev, nominalQ = TRUE) # rev

message("\n>> Saving error rate plots for forward and reverse trimmed reads on: ", 
        paste0("\n"), plot_dir, paste0("/fwd_error_rate_plot.png"),
        paste0("\n"), plot_dir, paste0("/rev_error_rate.png"), paste0("\n"))

# saving error rate plots for forward and reverse trimmed reads
ggsave(filename = "fwd_error_rate_plot.png", plot = fwd_plotErrors, path = plot_dir, width = 12, height = 7)
ggsave(filename = "rev_error_rate_plot.png", plot = rev_plotErrors, path = plot_dir, width = 12, height = 7)


#____Denoise_Unique_Sequences_______________________________________________________________________//

# The inference of Amplicon Sequence Variants (ASVs) was performed based on the previous estimated error rates. 
# In addition, the reads were dereplicated into unique sequences (100% identical) and singletons were removed.

message("\n>>> Denoising Sequences... \n")

message(">> Forward reads: ")
dadaFwd <- dada(derep = filtFastqFwdPath, err = errFwd, multithread = TRUE) # denoise fwd seqs
message("\n>> Reverse reads: ")
dadaRev <- dada(derep = filtFastqRevPath, err = errRev, multithread = TRUE) # denoise rev seqs

# inspect results
message("\n>> Inspect Forward reads: ")
dadaFwd[[1]]
message("\n>> Inspect Reverse reads: ")
dadaRev[[1]]


#____Merging_paired_reads_______________________________________________________________________//

message("\n>>> Merging Paired Reads... \n")

# Merge paired-end reads
mergePE <- mergePairs(dadaF = dadaFwd, derepF = filtFastqFwdPath, 
                      dadaR = dadaRev, derepR = filtFastqRevPath, verbose = TRUE)


#____Construct_ASV_table_______________________________________________________________________//

message("\n>>> Constructing the ASV table... \n")

# Make an ASV table
asvTbl <- makeSequenceTable(samples = mergePE) # tabulate ASVs

message("\n>> Saving the ASV table on: \n", 
        tables_dir, paste0("/asv_table.csv\n"))

# save table filterTrimReads
write.csv(asvTbl ,paste0(tables_dir,"/asv_table.csv"))


#____Remove_Chimeras___________________________________________________________________________//

message("\n>>> Removing Chimeras from the ASV table... \n")

# Remove chimeras from the ASV table
asvTblNoChim <- removeBimeraDenovo(unqs = asvTbl, method = "consensus",
                                   multithread = TRUE, verbose = TRUE) 

message("\n>> Saving the ASV table without Chimeras on: \n", 
        tables_dir, paste0("/asv_table_without_chimeras.csv\n"))

# save table filterTrimReads
write.csv(asvTblNoChim ,paste0(tables_dir,"/asv_table_without_chimeras.csv"))


#____Track_Reads_throught_the_pipeline___________________________________________________________________________//

# Summarize the no. of sequences kept in each pipeline step

# function that sums `sum(getUniques(x)` the no. of unique sequences `getUniques(x)`
getN <- function(x) sum(getUniques(x)) 

# build a matrix with all the sequences kept in each pipeline step
summaryTblSeq <- cbind(filterTrimReads, # initial reads and filtered/trimmed reads
                       sapply(dadaFwd, getN), sapply(dadaRev, getN), # denoised sequences 
                       sapply(mergePE, getN), # merged PE sequences
                       rowSums(asvTblNoChim)) # non-chimeric sequences

# rename the column and row names 
colnames(summaryTblSeq) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(summaryTblSeq) <- sampleNames

# percentage table
summaryTblSeqPerc <- apply(X = summaryTblSeq, MARGIN = 2, function(x) x / summaryTblSeq[,1] * 100)

message(">> Saving the table with absolute and percentage values that summarize the number of reads kept in each pipeline step on:")
message("\n> Absolute table: ", tables_dir, paste0("/summary_sequences_absolute.csv"))
message("> Percentage table: ", tables_dir, paste0("/summary_sequences_percentage.csv\n"))
 
# save table filterTrimReads
write.csv(summaryTblSeq, paste0(tables_dir,"/summary_sequences_absolute.csv"))
write.csv(summaryTblSeqPerc, paste0(tables_dir,"/summary_sequences_percentage.csv"))


#____Assign_Taxonomy______________________________________________________________________________________________________//

message("\n>>> Assignining Taxonomy against the non-redundant SILVA database v138... \n")

# naive Bayes classifier - # assign taxonomy against the SILVA database (version 138)
taxTbl <- assignTaxonomy(seqs = asvTblNoChim, 
                         refFasta = paste0(taxonomy_dir,"/database/silva_nr99_v138_train_set.fa.gz"), 
                         multithread = TRUE)

message("\n> Adding species into the previous assignment based on 100% match...\n")

# add species into the previous assignment based on 100% match
taxTbl <- addSpecies(taxtab = taxTbl, refFasta = paste0(taxonomy_dir,"/database/silva_species_assignment_v138.fa.gz"))


####----------------------------------_#####

# keep the trackability of ASVs 
# add a new column with the new ASV labels/ids to the taxonomy table
taxTbl2 <- cbind(taxTbl, "ASV" = paste0("ASV_", 1:nrow(taxTbl))) 

# substitute the DNA sequences in rownames by the new identifiers/tags/ids "ASV_nrSeq" in the taxonomy table
rownames(taxTbl2) <- taxTbl2[,8] 

# retrieve the DNA sequences 
uniquesToFasta(asvTblNoChim, paste0(tables_dir,"/asvFastaDNASequences.fasta"), ids = taxTbl2[,8])

# do the same for the ASV table (with the distribution)
asvTblNoChim2 <- asvTblNoChim # copy ASV table
colnames(asvTblNoChim2) <- taxTbl2[,8] # substitute column DNA sequences names by "ASV_nrSeq" in the ASV table
asvTblNoChim2 <- t(asvTblNoChim2) # transpose the ASV matrix table 
asvTblNoChim2 <- as.data.frame(asvTblNoChim2)
asvTblNoChim2[,"ASV_ID"] <- rownames(asvTblNoChim2)

message("\n> Saving ASV and Taxonomy tables in R format on:\n")
message("ASV table: ", paste0(taxonomy_dir,"/objects/asvTblNoChim.rds"))
message("Taxonomy table: ", paste0(taxonomy_dir,"/objects/taxTbl.rds\n"))

# save ASV and taxonomy tables in R format 
saveRDS(object = asvTblNoChim2, file = paste0(taxonomy_dir,"/objects/asvTblNoChim.rds")) # save ASV table
saveRDS(object = taxTbl2, file = paste0(taxonomy_dir,"/objects/taxTbl.rds")) # save taxonomy table

# put taxonomy in a compatible format to convert it latter to biom format
tax2biom <- function(taxTable) {
        
        lenNrow = nrow(taxTable)
        lenNcol = ncol(taxTable)
        id = c()
        fullTax = c()
        
        
        for (irow in 1:lenNrow) {
                rowID = as.character(rownames(taxTable)[irow])
                id = append(id, rowID)
                tax = ""
                
                for (icol in 1:lenNcol) {
                        colName <- colnames(taxTable)[icol]
                        colSign <- paste0(tolower(strtrim(colName,1)), "__")
                        
                        if (tax == "") {
                                tax <- paste0(colSign, taxTable[irow,icol])
                                #tax <- taxTable[irow,icol]
                        }
                        else {
                                tax <- paste(tax, paste0(colSign, taxTable[irow,icol]), sep = "; ")
                                #tax <- paste(tax, taxTable[irow,icol], sep = "; ")
                        }
                }
                fullTax <- append(fullTax, tax)
        }
        
        taxTable2Biom <- data.frame(id, fullTax)
        colnames(taxTable2Biom) <- c("ASV_ID", "taxonomy")
        
        return(taxTable2Biom)
}

taxTbl2 <- tax2biom(taxTbl2) 

# Join ASV and Taxonomy tables into one
## exclude the "ID" first column from "taxTbl2" because "asvTblNoChim2" has already this information
asvTaxTbl <- cbind(asvTblNoChim2, "taxonomy" = taxTbl2[,-1]) 

write.table(x = asvTaxTbl, file = paste0(tables_dir,"/asvTaxTbl.txt"), sep = "\t", row.names = FALSE, quote = FALSE) # save ASV-taxonomy tables
