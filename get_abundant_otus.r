otu.tab.taxonomy <- read.table("data/nash_data/otus_by_ALL.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE, stringsAsFactors=FALSE)

### PROCESS TABLE INTO NUMERIC AND EXTRACT TAXONOMY

otu.tab.taxonomy <- t(otu.tab.taxonomy)
taxonomy <- otu.tab.taxonomy[,ncol(otu.tab.taxonomy)]
otu.tab <- otu.tab.taxonomy[,c(1:(ncol(otu.tab.taxonomy)-1))]
otu.tab <- t(apply(otu.tab,1,as.numeric))
rownames(otu.tab) <- rownames(otu.tab.taxonomy)
colnames(otu.tab) <- colnames(otu.tab.taxonomy)[c(1:(ncol(otu.tab.taxonomy)-1))]

### CODE FOR FILTERING AND SELECTING SAMPLES

metagenomicNASH <- c("CL.166", "CL.169", "CL.139", "CL.173", "CL.177", "CL.160", "CL.165", "CL.119")
otu.metagenomicNASH <- otu.tab[,which(colnames(otu.tab)%in%metagenomicNASH)]
metagenomicHealthy <- c("HLD.100", "HLD.102", "HLD.111", "HLD.80", "HLD.85", "HLD.28", "HLD.112", "HLD.23")
otu.metagenomicHealthy <- otu.tab[,which(colnames(otu.tab)%in%metagenomicHealthy)]
otu.tab.metagenomic <- data.frame(nrow=(nrow(otu.tab)),ncol=(length(metagenomicNASH)+length(metagenomicHealthy)))
otu.tab.metagenomic <- cbind(otu.metagenomicNASH, otu.metagenomicHealthy)

### CODE FOR GETTING OTUS THAT ARE AT LEAST 0.2% ABUNDANT IN AT LEAST ONE SAMPLE

abundance.cutoff <- 0.002

otu.sums <- apply(otu.tab.metagenomic,2,sum)
one.percents <- otu.sums * abundance.cutoff

#determine which OTUs are > 1% abundance in any sample
which.otus.abundant <- rep(FALSE,length(rownames(otu.tab.metagenomic)))
for (i in 1:length(rownames(otu.tab.metagenomic))) {
  if (any(otu.tab.metagenomic[i,]>one.percents)) {
    which.otus.abundant[i] = TRUE
  }
}
abundant.otus <- rownames(otu.tab.metagenomic)[which.otus.abundant]

### CODE FOR OUTPUTTING FASTA FILE WITH ABUNDANT OTU SEQUENCES ONLY

fasta.in <- read.table("data/OTU_seed_seqs.fa", sep="\t",stringsAsFactors=FALSE)
fasta.in <- fasta.in[,1]
ids.string <- fasta.in[c(TRUE,FALSE)]
seqs <- fasta.in[c(FALSE,TRUE)]
ids <- gsub("^>lcl\\|[0-9]*\\|num\\|[0-9]*\\|OTU\\|","",ids.string)
fasta.out <- c(1:(length(abundant.otus)*2))
fasta.out[c(TRUE,FALSE)] <- ids.string[which(ids%in%abundant.otus)]
fasta.out[c(FALSE,TRUE)] <- seqs[which(ids%in%abundant.otus)]
write.table(fasta.out,file="data/OTU_seed_seqs_less_common_removed.fa",col.names=FALSE,row.names=FALSE,quote=FALSE)

### CODE FOR GETTING OTUS THAT ARE OVERALL AT LEAST 1% OR 0.1% ABUNDANT

otu.sum <- apply(otu.tab,1,sum)
total.count <- sum(otu.sum)

one.percent <- 0.01 * total.count
abundant.otus <- rownames(otu.tab)[which(otu.sum > one.percent)]

point.one.percent <- 0.001 * total.count
less.abundant.otus <- rownames(otu.tab)[which(otu.sum > point.one.percent)]

abundant.sum <- sum(otu.sum[which(rownames(otu.tab)%in%abundant.otus)])
abundant.sum/total.count
# [1] 0.5860075

less.abundant.sum <- sum(otu.sum[which(rownames(otu.tab)%in%less.abundant.otus)])
less.abundant.sum/total.count
# [1] 0.954467

length(less.abundant.otus)
# [1] 140

exclude.otus <- rownames(otu.tab)[which(!rownames(otu.tab)%in%less.abundant.otus)]

# manually removed exclude OTUs from OTU seed sequence list (OTU_seed_seqs.fa)
