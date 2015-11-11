#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

input_file = args[1]
output_histogram_file = args[2]
output_sequence_length_over_cutoff_file = args[3]

### for running within R
# input_file <- "data/orfLengths.txt"
# output_histogram_file <- "data/orf_length_histogram.pdf"
# output_sequence_length_over_cutoff_file <- "data/too_long_orfs.txt"

print(paste("reading file",input_file,"(this may take a while if your file is large)"))
data <- read.table(input_file,sep="\t",header=FALSE)
colnames(data) <- c("filename","seqname","seqlength")

print("sorting sequences by length")
# sort ORFs by ascending sequence length
data <- data[order(data$seqlength),]

low_cutoff <- 60
high_cutoff <- 3000

print(paste("printing sequence length histogram to PDF",output_histogram_file))
pdf(output_histogram_file)
hist(log2(data$seqlength),main="frequency of lengths of predicted orfs",xlab="log base 2 of sequence length")
hist(data$seqlength,breaks=c(seq(0,high_cutoff,100),max(data$seqlength)),freq=FALSE,main="density of lengths of predicted orfs",xlab=paste("sequence length with cutoff of",high_cutoff),xlim=c(0,high_cutoff))
hist(data$seqlength,breaks=c(seq(0,high_cutoff*2,100),max(data$seqlength)),freq=FALSE,main="density of lengths of predicted orfs",xlab=paste("sequence length with cutoff of",(high_cutoff*2)),xlim=c(0,high_cutoff*2))
dev.off()

print(paste("writing sequence length table for sequences over cutoff to file",output_sequence_length_over_cutoff_file))
high_outliers <- data[which(data$seqlength>high_cutoff),]
low_outliers <- data[which(data$seqlength<low_cutoff),]
outliers <- rbind(high_outliers,low_outliers)
write.table(outliers,file=output_sequence_length_over_cutoff_file,sep="\t",quote=FALSE,row.names=FALSE)

print(paste(nrow(outliers),"predicted coding sequences were dropped out of",nrow(data)," for a drop rate of",(nrow(outliers)/nrow(data)),"percent."))
