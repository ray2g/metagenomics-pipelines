## João Raimundo @BioData.pt ----------------------///

#______Import_Packages________________________//

message('\nDADA2 Pipeline\nJoão Raimundo @BioData.pt\n')
message('\n>>> Importing Packages... \n')

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

message("\n>>> Reading data... \n") 
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

message("\n>> Saving read quality profile plots on:\n")
message(plot_dir, paste0("/fwd_qualityProfile.png"))
message(plot_dir, paste0("rev_qualityProfile.png"), paste0("\n"))

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

message(">> Saving trimmed reads quality profile plots on:\n") 
message(plot_dir, paste0("/fwd_qualityProfile_trimmed.png"))
message(plot_dir, paste0("/rev_qualityProfile_trimmed.png"), paste0("\n"))


# saving trimmed reads quality profiles plots
ggsave(filename = "fwd_qualityProfile_trimmed.png", plot = fwd_trimmed_qualityProfile, path = plot_dir, width = 12, height = 8)
ggsave(filename = "rev_qualityProfile_trimmed.png", plot = rev_trimmed_qualityProfile , path = plot_dir, width = 12, height = 8)

message(">> Saving the trimmed summary table on:\n") 
message(tables_dir, paste0("/number_sequences_before_and_after_trimming.tsv\n"))

# save table filterTrimReads
write.table(x = filterTrimReads, file = paste0(tables_dir,"/number_sequences_before_and_after_trimming.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) 

#____Learning_the_Error_Rates_________________________________________________________________________________________//

message(">>> Learning the Error Rates... \n")

# model/learn the error rates for the filtered fastq files
message("> Forward reads:")
errFwd <- learnErrors(fls = filtFastqFwdPath, multithread = TRUE) # fwd

message("\n> Reversed reads:")
errRev <- learnErrors(fls = filtFastqRevPath, multithread = TRUE) # rev

# Plot errors 
fwd_plotErrors <- plotErrors(dq = errFwd, nominalQ = TRUE) # fwd
rev_plotErrors <- plotErrors(dq = errRev, nominalQ = TRUE) # rev

message("\n>> Saving error rate plots for forward and reverse trimmed reads on:")
message(paste0("\n"), plot_dir, paste0("/fwd_error_rate_plot.png"))
message(plot_dir, paste0("/rev_error_rate.png"), paste0("\n"))

# saving error rate plots for forward and reverse trimmed reads
ggsave(filename = "fwd_error_rate_plot.png", plot = fwd_plotErrors, path = plot_dir, width = 12, height = 7)
ggsave(filename = "rev_error_rate_plot.png", plot = rev_plotErrors, path = plot_dir, width = 12, height = 7)


#____Denoise_Unique_Sequences_______________________________________________________________________//

# The inference of Amplicon Sequence Variants (ASVs) was performed based on the previous estimated error rates. 
# In addition, the reads were dereplicated into unique sequences (100% identical) and singletons were removed.

message("\n>>> Denoising Sequences... \n")

message("> Forward reads: ")
dadaFwd <- dada(derep = filtFastqFwdPath, err = errFwd, multithread = TRUE) # denoise fwd seqs
message("\n> Reverse reads: ")
dadaRev <- dada(derep = filtFastqRevPath, err = errRev, multithread = TRUE) # denoise rev seqs

# inspect results
message("\n>> Inspect Forward reads: ") # fwd
dadaFwd[[1]]

message("\n>> Inspect Reverse reads: ") #rev 
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

message("\n>> Saving the ASV table on: \n", tables_dir, paste0("/asv_table.tsv\n"))

# save table filterTrimReads
write.table(x = asvTbl, file = paste0(tables_dir,"/asv_table.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) 

#____Remove_Chimeras___________________________________________________________________________//

message("\n>>> Removing Chimeras from the ASV table... \n")

# Remove chimeras from the ASV table
asvTblNoChim <- removeBimeraDenovo(unqs = asvTbl, method = "consensus",
                                   multithread = TRUE, verbose = TRUE) 

message("\n>> Saving the ASV table without Chimeras on: \n", 
        tables_dir, paste0("/asv_table_without_chimeras.tsv\n"))

# save table filterTrimReads
write.table(x = asvTblNoChim, file = paste0(tables_dir,"/asv_table_without_chimeras.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) 


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
message("\n> Absolute table: ", tables_dir, paste0("/summary_sequences_absolute.tsv"))
message("> Percentage table: ", tables_dir, paste0("/summary_sequences_percentage.tsv\n"))
 
# save table filterTrimReads
write.table(x = summaryTblSeq, file = paste0(tables_dir,"/summary_sequences_absolute.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) 

write.table(x = summaryTblSeqPerc, file = paste0(tables_dir,"/summary_sequences_percentage.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) 


#____Assign_Taxonomy______________________________________________________________________________________________________//

message(">>> Assignining Taxonomy against the non-redundant SILVA database v138... \n")

# naive Bayes classifier - # assign taxonomy against the SILVA database (version 138)
taxTbl <- assignTaxonomy(seqs = asvTblNoChim, 
                         refFasta = paste0(taxonomy_dir,"/database/silva_nr99_v138_train_set.fa.gz"), 
                         multithread = TRUE)

message("\n>>> Adding species into the previous assignment based on 100% match...\n")

# add species into the previous assignment based on 100% match
taxTbl <- addSpecies(taxtab = taxTbl, refFasta = paste0(taxonomy_dir,"/database/silva_species_assignment_v138.fa.gz"))


#___keep_the_trackability_of_ASVs_______________________________________________________________________________________/

# add a new column with the new ASV labels/ids to the taxonomy table
taxTbl2 <- cbind("ASV" = paste0("ASV_", 1:nrow(taxTbl)),taxTbl)

# substitute the DNA sequences in rownames by the new identifiers/tags/ids "ASV_nrSeq" in the taxonomy table
rownames(taxTbl2) <- taxTbl2[,1] 

message(">> Saving ASV DNA Sequences on:\n")
message(paste0(path ,"asvFastaDNASequences.fasta"))

# retrieve the DNA sequences 
uniquesToFasta(asvTblNoChim, paste0(path ,"asvFastaDNASequences.fasta"), ids = taxTbl2[,1])

# do the same for the ASV table (with the distribution)
asvTblNoChim2 <- asvTblNoChim # copy ASV table
colnames(asvTblNoChim2) <- taxTbl2[,1] # substitute column DNA sequences names by "ASV_nrSeq" in the ASV table
asvTblNoChim2 <- t(asvTblNoChim2) # transpose the ASV matrix table 
asvTblNoChim2 <- as.data.frame(asvTblNoChim2)
ASV_ID <- rownames(asvTblNoChim2)
asvTblNoChim2 <- cbind(ASV_ID, asvTblNoChim2)

message("\n>> Saving ASV and Taxonomy tables in tsv format on:\n")
message("> ASV table: ", paste0(tables_dir,"/taxTbl.tsv"))
message("> Taxonomy table: ", paste0(tables_dir,"/taxTbl.tsv"))

# save ASV and taxonomy tables in TSV format
write.table(x = asvTblNoChim2, file = paste0(tables_dir,"/asvTblNoChim.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) # save ASV table

write.table(x = taxTbl2, file = paste0(tables_dir,"/taxTbl.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE) # save taxonomy table 

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
        colnames(taxTable2Biom) <- c("ASV_ID", "Taxonomy")
        
        return(taxTable2Biom)
}

taxTbl2 <- tax2biom(taxTbl2) 

message("\n\n>>> Joining ASV and Taxonomy tables into one...\n")

# Join ASV and Taxonomy tables into one - exclude the "ID" first column from "taxTbl2" 
# because "asvTblNoChim2" has already this information
asvTaxTbl <- cbind(asvTblNoChim2, "Taxonomy" = taxTbl2[,-1]) 

message(">> Saving the joined ASV-Taxonomy table in TSV format on:\n")
message(paste0(tables_dir,"/asvTaxTbl.tsv"))

# save ASV-taxonomy table
write.table(x = asvTaxTbl, file = paste0(tables_dir,"/asvTaxTbl.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("\n\n>>> Converting ASV-Taxonomy table into BIOM format...\n")

# convert asvTaxTbl table into biom - function
convertTab2Biom <- function(inFile, outFile) {
  
  if (system("command -v biom", ignore.stdout = TRUE, ignore.stderr = TRUE) !=0)  {
    
    stop("biom program is not installed or it is not accessible!\n  Exiting...")
    
  }
  
  system(paste("biom convert", "-i", inFile, "-o", outFile, "--to-hdf5", 
               '--table-type="OTU table"', "--process-obs-metadata taxonomy"))
  
}

# convert ASV table with taxonomy in tab-delimited format into biom format  
convertTab2Biom(inFile = paste0(tables_dir,"/asvTaxTbl.tsv"), outFile = paste0(tables_dir,"/asvTaxTbl.biom"))

message(">> Saving the joined ASV-Taxonomy table in BIOM format on:\n")
message(paste0(tables_dir,"/asvTaxTbl.biom"))
                           
                           
