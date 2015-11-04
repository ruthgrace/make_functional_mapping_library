#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

input_file = args[1]
output_sequence_length_file = args[2]
output_histogram_file = args[3]
output_sequence_length_over_cutoff_file = args[4]

print(paste("reading file",input_file,"(this may take a while if your file is large)"))
data <- readLines(input_file)
print("creating sequence length table")
df <- data.frame(matrix(ncol=3,nrow=(length(data)/2)))
row <- 1
filename <- ""
seqname <- ""
colnames(df) <- c("filename","sequence_name","sequence_length")
print("adding sequences to sequence length table")
for (line in data) {
  if (grepl("\\.fa$",line)) {
    filename <- line
  }
  else if (grepl("^>",line)) {
    seqname <- line
  }
  # skip if line is empty
  else if (line != ""){
    df[row,] <- c(filename, seqname, line)
    row <- row + 1
  }
}

print("sorting sequences by length")
# sort ORFs by ascending sequence length
df <- df[order(df$sequence_length),]

print(paste("writing sequence length table to file",output_sequence_length_file))
write.table(df,file=output_sequence_length_file,sep="\t",quote=FALSE)

cutoff <- 5000

print(paste("printing sequence length histogram to PDF",output_histogram_file))
pdf(output_histogram_file)

histogram(df$sequence_length,main="frequency of lengths of predicted orfs")

dev.off()

print(paste("writing sequence length table for sequences over cutoff to file",output_sequence_length_over_cutoff_file))
write.table(df[which(df$sequence_length>cutoff),],file=output_sequence_length_over_cutoff_file,sep="\t")