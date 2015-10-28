#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

print(args)

input_file = args[1]
output_file = args[2]

print(paste("sorting sequences by closeness to median length in file",input_file,"output in file",output_file))

fasta <- read.table(input_file,sep="\t",header=FALSE)
fasta <- fasta[,1]
ids <- fasta[c(TRUE,FALSE)]
seqs <- fasta[c(FALSE, TRUE)]
seqs <- as.character(seqs)
seq_lengths <- nchar(seqs)
median <- median(seq_lengths)
median_dist <- abs(seq_lengths - median)
sorted_order <- order(median_dist)
ids <- ids[sorted_order]
seqs <- seqs[sorted_order]

sorted_fasta <- data.frame(matrix(ncol=1,nrow=(length(ids) + length(seqs))))

sorted_fasta[c(TRUE,FALSE),1] <- as.character(ids)
sorted_fasta[c(FALSE,TRUE),1] <- as.character(seqs)

write.table(sorted_fasta,file=output_file,row.names=FALSE,col.names=FALSE,quote=FALSE)
