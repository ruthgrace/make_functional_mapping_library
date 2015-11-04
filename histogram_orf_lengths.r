print("reading file (this may take a while if your file is large)")
data <- readLines("count_orf_lengths_nohup.out")
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

print("writing sequence length table to file")
write.table(df,file="data/nucleotides_per_count_per_sequence_per_file.txt",sep="\t",quote=FALSE)

cutoff <- 5000

print("printing sequence length histogram to PDF")
pdf("data/orf_length_histogram.pdf")

histogram(df$sequence_length,main="frequency of lengths of predicted orfs")

dev.off()

print("writing sequence length table for sequences over cutoff to file")
write.table(df[which(df$sequence_length>cutoff),],file="data/too_long_orfs.txt",sep="\t")