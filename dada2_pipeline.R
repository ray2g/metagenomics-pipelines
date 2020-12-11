## Joao Raimundo @BioData.pt ----------------------///


#______Import_Packages________________________//

message('\n Importing Packages... \n')

library("dada2") 
library("ggplot2")


#______Directories______________________________________________________________________________________________//

setwd("~/Desktop/dada2_test/") # set the working directory ## CHANGE ME

path <- "~/Desktop/dada2_test/data/" # set the data directory ## CHANGE ME  

plot_dir <- "~/Desktop/dada2_test/plots" # set the directory to save the plots ## CHANGE ME

fastqPath <- "~/Desktop/dada2_test/data/fastq" # set the directory to save processed fastq files ## CHANGE ME

tables_dir <- "~/Desktop/dada2_test/tables" # set the directory to save the tables ## CHANGE ME

message("\n Reading data... \n") 
list.files(path)


#____Match_Forward_and_Reverse_fastq_files________________________________________________________________________//

#set the relative path to to each forward and reverse fastq files directory
fastqFwdPath <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fastqRevPath <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sampleNames <- sapply(strsplit(basename(fastqFwdPath), "_"), `[`, 1)


#____Inspect_read_quality_profiles________________________________________________________________________________//

message("\n Inspecting Read Quality Profiles... \n")

fwd_qualityProfile <- plotQualityProfile(fl = fastqFwdPath) # forward 
rev_qualityProfile <- plotQualityProfile(fl = fastqRevPath) # reverse

message("\n Saving read quality profile plots on: ", paste0("\n"), plot_dir, paste0("/fwd_qualityProfile.png"), 
                                                     paste0("\n"), plot_dir, paste0("rev_qualityProfile.png"), paste0("\n"))

# saving quality profiles plots
ggsave(filename = "fwd_qualityProfile.png", plot = fwd_qualityProfile, path = plot_dir, width = 12, height = 8)
ggsave(filename = "rev_qualityProfile.png", plot = rev_qualityProfile , path = plot_dir, width = 12, height = 8)


#___Filter_and_trim_reads______________________________________________________________________________________________________________//

message("\n Filter and Trimming the Reads... \n")
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

message(paste0("\n The filtered reads were saved on: \n"),fastqPath,paste0("/filtered/\n"))

message("Number of sequences kept after filtering/trimming in relation to the initial number of reads: \n")
filterTrimReads

message("\n")

message("\n Inspecting Trimmed Read Quality Profiles... \n")

# forward quality plots after filtering and trimming
fwd_trimmed_qualityProfile <- plotQualityProfile(fl = filtFastqFwdPath)

# reverse quality plots after filtering and trimming
rev_trimmed_qualityProfile <- plotQualityProfile(fl = filtFastqRevPath)

message("\n Saving trimmed reads quality profile plots on: ", 
        paste0("\n"), plot_dir, paste0("/fwd_qualityProfile_trimmed.png"), 
        paste0("\n"), plot_dir, paste0("/rev_qualityProfile_trimmed.png"), paste0("\n"))

# saving trimmed reads quality profiles plots
ggsave(filename = "fwd_qualityProfile_trimmed.png", plot = fwd_trimmed_qualityProfile, path = plot_dir, width = 12, height = 8)
ggsave(filename = "rev_qualityProfile_trimmed.png", plot = rev_trimmed_qualityProfile , path = plot_dir, width = 12, height = 8)

message("\n Saving the trimmed summary table on: \n", 
        tables_dir, paste0("/number_sequences_before_and_after_trimming.csv\n"))

# save table filterTrimReads
write.csv(filterTrimReads,paste0(tables_dir,"/number_sequences_before_and_after_trimming.csv"))


#____Learning_the_Error_Rates_______________________________________________________________________________________________________________//

message("\n Learning the Error Rates... \n")

# model/learn the error rates for the filtered fastq files
errFwd <- learnErrors(fls = filtFastqFwdPath, multithread = TRUE) # fwd
errRev <- learnErrors(fls = filtFastqRevPath, multithread = TRUE) # rev

# Plot errors 
fwd_plotErrors <- plotErrors(dq = errFwd, nominalQ = TRUE) # fwd
rev_plotErrors <- plotErrors(dq = errRev, nominalQ = TRUE) # rev

message("\n Saving error rate plots for forward and reverse trimmed reads on: ", 
        paste0("\n"), plot_dir, paste0("/fwd_error_rate_plot.png"),
        paste0("\n"), plot_dir, paste0("/rev_error_rate.png"), paste0("\n"))

# saving error rate plots for forward and reverse trimmed reads
ggsave(filename = "fwd_error_rate_plot.png", plot = fwd_plotErrors, path = plot_dir, width = 12, height = 7)
ggsave(filename = "rev_error_rate_plot.png", plot = rev_plotErrors, path = plot_dir, width = 12, height = 7)


