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

I wanted all the genomes that matched 98% or higher for each OTU, or the best match if there were no matches better than 98%. The script extracts the gi numbers for such matches, and then crawls the NCBI nuccore website to find the taxon ID. The taxon ID is found in ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt and the corresponding FTP link is used to download the genome. For each species found by this method, the genomes for 10 random strains are downloaded (or all of the strains if there are less than 10). The idea with having 10 representatives is to increase the coverage of the library.

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

As a sanity check, you may want to throw a few of your ORFs into a protein translator, and blast them to make sure they're real protein sequences (no frame shift mutations or anything).

## Checking ORF length

This script concatenates all the fasta files in each genus folder in `data/orfs/` into one fasta file per genus.

```bash
nohup ./concatenate_orfs_by_genus.sh > concatenate_orfs_nohup.out 2>&1&
```

To account for incorrectly annotated ORFs, run the following program to collect all ORF lengths into the tab separated file `data/orfLengths.txt`.

```bash
nohup ./count_orf_lengths.sh > data/orfLengths.txt 2>&1&
```

Plot the lengths on a histogram by running the `histogram_orf_lengths.r` script. The histogram is output into `data/orf_length_histogram.pdf`.

```R
nohup Rscript histogram_orf_lengths.r "data/orfLengths.txt" "data/orf_length_histogram.pdf" "data/too_long_orfs.txt" > histogram_orf_lengths_nohup.out 2>&1&
```

Look at the histogram and decide what coding sequence length cutoff is reasonable for differentiating real ORFs and too long artifacts. Change this line in the `histogram_orf_lengths.r` script to your cutoff of choice, run the script again, and the sequences that are too long will be output into `data/too_long_orfs.txt`. The default cutoff is 3000:

```R
low_cutoff <- 60
high_cutoff <- 3000
```

Remove or correct the outlier sequences before clustering. To automatically remove all sequences that don't make the cutoff, run:

```bash
nohup ./remove_orfs.pl "data/too_long_orfs.txt" > remove_orfs_nohup.out 2>&1&
```

## Clustering at 99% by genus for BLAST

The orfs for each genus are then sorted for deviation from median length, and clustered at 99% identity using CD-HIT. CD-HIT is fast and doesn't use too much memory because it doesn't do real centroid picking (it picks the first sequences as seeds), and the best centroids for 100% identity clustering are the orfs with median length. You can download CD-HIT [here](https://code.google.com/p/cdhit/downloads/detail?name=cd-hit-v4.6.1-2012-08-27.tgz) and install instructions are [here](http://weizhong-lab.ucsd.edu/cd-hit/wiki/doku.php?id=cd-hit_user_guide). I followed the instructions for a multithreaded install.

```bash
nohup ./cluster_orfs_by_genus_99_multithreaded.sh > cluster_orfs_by_genus_99_nohup.out 2>&1&
```

The resulting sequences are what we will be blasting to the SEED annotation database. Sequences will be mapped using Bowtie2, which won't work properly if there are spaces in the sequence identifiers.

To remove the space between the ">" and the unique number for each identifier, run:

```
sed -i.backup '/^>/ s/^> />/' all_genus_orfs_clustered_at_99_unique.fa 
```

To replace the rest of the spaces in the sequence identifiers with underscores, run:

```
sed -i.backup '/^>/ s/ /_/g' all_genus_orfs_clustered_at_99_unique.fa 
```

## Assigning function

Add the genus name to each sequence ID just in case:

```bash
nohup ./add_genus_name_to_seq_ids.sh > add_genus_name_to_seq_ids_nohup.out 2>&1&
```

Concatenate all the 99% clustered per genus sequences into a single file:

```bash
nohup cat data/orfs/*/*_multithreaded_99_cd_hit_genus_id.txt > data/orfs/all_genus_orfs_clustered_at_99.fa 2>&1&
```

Prepend unique number to beginning of sequence identifier, to make all identifiers unique, just in case:

```bash
nohup ./uniqueify_seq_id.pl data/orfs/all_genus_orfs_clustered_at_99.fa data/orfs/all_genus_orfs_clustered_at_99_unique.fa > uniqueify_seq_id_nohup.out 2>&1&
```

The sequences then have to be translated into protein. I'm using Biopython. Installation instructions are [here](http://biopython.org/DIST/docs/install/Installation.html).

```bash
nohup python translate_sequences.py data/orfs/all_genus_orfs_clustered_at_99_unique.fa data/orfs/all_genus_orfs_clustered_at_99_unique_protein.fa
```

Check for spurious stop codons to filter out incorrectly annotated proteins (denoted asterisks and dashes). If you don't have any, this line won't output anything.

```bash
grep "^[^>].*[*-].*" data/orfs/all_genus_orfs_clustered_at_99_unique_protein.fa
```

If they exist you can remove them like so:

```bash
nohup ./remove_incorrect_proteins.pl data/orfs/all_genus_orfs_clustered_at_99_unique_protein.fa data/orfs/all_genus_orfs_clustered_at_99_unique_protein_validated.fa data/orfs/removed_protein_sequences.fa > remove_incorrect_proteins_nohup.out 2>&1&
```

Run BLAST with the SEED database.

Database details from Jean: "SEED databse downloaded June 2013. I also added in missing fig.pegs from the old SEED database (2010) that we had used for the Microbiome paper (Macklaim et al. 2013)"

The command and run time for an example BLAST:

This was for 413986 sequences

```
nohup blastp -db db_fastas.complex.faa -query leftover_refseqs_blast_seed.faa -out leftover_refseqs_blast_seed.faa.out -outfmt 6 -evalue 1e-3 -num_alignments 10 -num_threads 4 > leftover_nohup.out 2>&1&
```

Jul 5, 2013 9:00am 25% done
Jul 7, 2013 7:30pm 61%
July 8, 2103 10:34am 70%
Jul 9, 2013 10:41am 85%

## Mapping

Mapping is done with the script in the mapping folder. Be sure to change the variables `WORKING_DIR`, `REFSEQS`, `OUTPUTFOLDER`, `IDX`, and `BIN` appropriately. Be careful to save your nohup file because this is your mapping index.

```
nohup mapping_scripts/mapping.sh > mapping_nohup.out 2>&1&
```

You may find it advantageous to split up the first and second half of the script, especially if the second half isn't working. I debugged the second half of the script for my specific use case in `mapping_count_aggregation_only.sh`.

#Things I tried that didn't work

Originally the plan was to cluster at 100% by genus, then 95% across all genus, and then blast to assign function. However, the 100% clustering resulted in about 4.5 million sequences in my case (after filtering out sequences longer than 3000). On our computer with 64GB of RAM and 16 threads, it would probably have taken a month to run the 95% clustering, which was too slow.

We then tried clustering by genus at 90%, which was too slow and only seemed to reduce the number of sequences ~87%. Similarly clustering by 95% per genus is slow - it took 15 minutes to do one genus (i have a couple hundred)

I also intended to map to the 100% within genus clustering and BLAST to the 99% within genus clustering, and amalgamate the mapped counts to the 99% clustering, but I actually just mapped to 99% within genus clustering and no amalgamation was neccessary. I am getting > 50% mapped reads on average. It's worth exploring what difference mapping to my 100% within genus clustering would have made.

