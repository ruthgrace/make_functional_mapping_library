# Collect Reference Genomes

This is Ruth's case-specific workflow for finding genomes to create a functional mapping library for a metagenomic study. The steps here include that you have R and Perl installed.

## Find out which OTUs are >= 0.2% abundance

In this experiment, the sequencing depth will have the power to detect a 2 fold change up or down in bacteria that are 0.2% abundant in a sample. The code for extracting the relevant OTU sequences is in get_abundant_otus.r

The OTU table is read in this line

```R
otu.tab.taxonomy <- read.table("data/nash_data/otus_by_ALL.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE, stringsAsFactors=FALSE)
```
and the OTU sequences are read in this line

```R
fasta.in <- read.table("data/OTU_seed_seqs.fa", sep="\t",stringsAsFactors=FALSE)
```

The different sections of code in get_abundant_otus.r do different things, so cut and paste what you need.

## BLAST OTU seed sequences with Human Microbiome Project data set

The gut HMP reference genomes were downloaded from http://hmpdacc.org/HMRGD/healthy/

I turned the .fa file into a local BLAST database
```
makeblastdb -in all_seqs.fa -dbtype 'nucl' -out all_seqs.fa
```

Run local BLAST with the OTU seed sequences as a query and the HMP reference as a database:

```bash
nohup blastn -db hmp_genomes/all_seqs.fa -query data/OTU_seed_seqs_less_common_removed.fa -out output/blast.out -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -evalue 1e-3 -num_alignments 10 -num_threads 4 > blast_nohup.out 2>&1&
```

Add this header to blast.out.
```
qseqid  sseqid  pident  length  mismatch  gapopen qstart  qend  sstart  send  evalue  bitscore  qlen  slen
```

Output list of sequences that don’t have a > 97% identity match using get_otus_not_in_hmp.r.

Note: I have decided that I will actually skip this step and just bring in all the HMP genomes later along with all the matching genomes I find through BLAST.

## Find matches in the NCBI complete and draft bacterial genomes for OTUs

I ran into a bug using the NCBI webtool, and had to search once through the wgs database, and once with Complete Genomes to get both the draft and the complete genomes.

Blast against draft and complete bacterial genomes with NCBI: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&BLAST_SPEC=MicrobialGenomes

I put my BLAST output (Downloaded in CSV format) in data/complete-genomes-Alignment-HitTable.csv and data/wgs-Alignment-HitTable.csv.

## Retrieve the genome for the matches

This is done in the script extract_genomes_not_in_hmp.py. It requires you to make a folder named genome inside the folder named data, and uses the blast outputs generated in the previous section.

I only wanted genomes that were the best match for each OTU. If there was a match > 99% percent identity, then I wanted a second genome that had a percent identity just over 97%. The script extracts the gi numbers for such matches, and then crawls the NCBI nuccore website to find the taxon ID. The taxon ID is found in ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt and the corresponding FTP link is used to download the genome.

The script can be run on the command line:

```bash
nohup python extract_genomes.py > extract_genomes_nohup.out 2>&1&
```

The download genomes are compressed. To uncompress, navigate to the folder in terminal and run:
```bash
gunzip *.gz
```

The end result is that all the genomes I need for my functional mapping library are in data/genomes/. For each genome, I've downloaded a feature table and a .fna file, so that the coding sequences can be extracted in the next step. Some of the sequences do not have a feature table, and in this case I've downloaded the .gff instead.

//TODO: Make sure that genus match even if surrounded by square brackets (e.g. "[Clostridium]" and "Clostridium" are the same genus)

## Extract coordinates of coding sequences from genomes

ORFS for .gff files need to be extracted using Glimmer. `extract_orfs.sh` is a bash script that runs the first steps of the iterated Glimmer script on all genomes for which we were only able to colleted a .gff file and not a feature_table.txt. The coordinates are output in the data/glimmer_output folder. Glimmer can be downloaded from the [John's Hopkins Center for Computational Biology website](https://ccb.jhu.edu/software/glimmer/).

```bash
nohup ./extract_orf_coordinates_from_gff.sh > extract_orfs_nohup.out 2>&1&
```

Note that after I ran this, a couple of output files were empty, presumably because no ORFs were found.

## Output coding sequences into fasta files

Note that after this step, genus names in brackets (such as "[Clostridium]") will be grouped with species name without brackets. Brackets denote uncertainty about the genus designation (the isolate has been named by sequence identity rather than biochemical evidence).

This bash script uses the `get_orf_sequences_from_feature_table.pl` Perl script to write the ORF sequences denoted by the feature tables into a fasta file, into a folder named the same as the genus in data/orfs/.

```bash
nohup ./output_feature_table_orf_sequences.sh > feature_table_orfs_output_nohup.out 2>&1&
```

This bash script uses the `get_orf_sequences_from_glimmer_output.pl` Perl script to write the ORF sequences denoted by the .coord files output by Glimmer into a fasta file, into a folder named the same as the genus in data/orfs/.

```bash
nohup ./output_glimmer_coord_orf_sequences.sh > glimmer_coord_orfs_output_nohup.out 2>&1&
```

As a sanity check, you may want to throw your ORFs into a protein translator, and blast them to make sure they're real protein sequences.

## Clustering at 100% by genus

This script concatenates all the fasta files in each genus folder in `data/orfs/` into one fasta file per genus.

```bash
nohup ./concatenate_orfs_by_genus.sh > concatenate_orfs_nohup.out 2>&1&
```

The orfs for each genus are then sorted for deviation from median length, and clustered at 100% identity using CD-HIT. CD-HIT is fast and doesn't use too much memory because it doesn't do real centroid picking (it picks the first sequences as seeds), and the best centroids for 100% identity clustering are the orfs with median length. You can download CD-HIT [here](https://code.google.com/p/cdhit/downloads/detail?name=cd-hit-v4.6.1-2012-08-27.tgz) and install instructions are [here](http://weizhong-lab.ucsd.edu/cd-hit/wiki/doku.php?id=cd-hit_user_guide).

```bash
nohup ./cluster_orfs_by_genus.sh > cluster_orfs_by_genus_nohup.out 2>&1&
```

## Clustering at 95% between genus

Concatenate all the 100% clustered per genus sequences into a single file:

```
cat data/orfs/*/*_cd_hit.txt > data/orfs/all_genus_orfs_clustered_at_100.fa
```

Cluster by 95% identity across genus using CD-HIT

```
nonhup cd-hit-v4.6.1-2012-08-27/cd-hit -i data/orfs/all_genus_orfs_clustered_at_100.fa -o data/orfs/all_orfs_clustered_at_95.txt -c 0.95 -n 5 > cluster_all_orfs_at_95_nohup.out 2>&1&
```
