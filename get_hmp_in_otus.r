# FIX THIS SCRIPT TO EXTRACT WHICH HMP REFERENCE GENOMES MATCH TO OTUS

# i guess what you really want is binomial species name

# add this header to blast.out:
# qseqid  sseqid  pident  length  mismatch  gapopen qstart  qend  sstart  send  evalue  bitscore  qlen  slen

data <- read.table("output/blast.out",header=TRUE)
otus.in.hmp <- data[which(data$pident >= 97),]
otus.not.in.hmp <- unique(data$qseqid[which(!as.character(data$qseqid)%in%as.character(otus.in.hmp$qseqid))])

fasta <- read.table("data/OTU_seed_seqs_less_common_removed.fa",sep="\t",header=FALSE)
fasta <- fasta[,1]
ids <- fasta[c(TRUE,FALSE)]
ids <- gsub(">","",as.character(ids))
seqs <- fasta[c(FALSE, TRUE)]

ids.not.in.hmp <- ids[which(as.character(ids)%in%as.character(otus.not.in.hmp))]
seqs.not.in.hmp <- seqs[which(as.character(ids)%in%as.character(otus.not.in.hmp))]

fasta.not.in.hmp <- data.frame(matrix(ncol=1,nrow=(length(ids.not.in.hmp) + length(seqs.not.in.hmp))))

ids.not.in.hmp <- gsub("^",">",ids.not.in.hmp)

fasta.not.in.hmp[c(TRUE,FALSE),1] <- as.character(ids.not.in.hmp)
fasta.not.in.hmp[c(FALSE,TRUE),1] <- as.character(seqs.not.in.hmp)

write.table(fasta.not.in.hmp,file="data/OTUs_not_in_HMP.fa",row.names=FALSE,col.names=FALSE,quote=FALSE)
