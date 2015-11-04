#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

input_file = args[1]
output_histogram_file = args[2]
output_sequence_length_over_cutoff_file = args[3]

print(paste("reading file",input_file,"(this may take a while if your file is large)"))
data <- read.table(input_file,sep="\t",header=FALSE)
colnames(data) <- c("filename","seqname","seqlength")

print("sorting sequences by length")
# sort ORFs by ascending sequence length
data <- data[order(data$seqlength),]

cutoff <- 3000

print(paste("printing sequence length histogram to PDF",output_histogram_file))
pdf(output_histogram_file)
hist(log2(data$seqlength),main="frequency of lengths of predicted orfs",xlab="log base 2 of sequence length")
dev.off()

print(paste("writing sequence length table for sequences over cutoff to file",output_sequence_length_over_cutoff_file))
outliers <- data[which(data$seqlength>cutoff),]
write.table(outliers,file=output_sequence_length_over_cutoff_file,sep="\t",quote=FALSE,row.names=FALSE)

print(paste(nrow(outliers),"predicted coding sequences were dropped out of",nrow(data)," for a drop rate of",(nrow(outliers)/nrow(data)),"percent."))
