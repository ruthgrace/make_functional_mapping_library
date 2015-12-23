#!/usr/bin/env bash

# June 14, 2013 - JM
# run: nohup ./bin/map_bac.sh > bowtie2_log_`date +"%Y-%m-%d_%H%M%S"`.txt 2>&1&

#------------------------------------------------------------------------------------------------------

# Let's start mapping to bacteria db
# These seqs are clustered at 90% id by uclust_smallmem. The seqs were ordered by abs deviation from median gene length

# Samples are in: /Volumes/rhamnosus/twntyfr/samples/
#	named: 001A.nonhuman.fastq etc.

# All mapping output will be put to WORKING_DIR as SAMPLE_NAME_map/
#WORKING_DIR=/Volumes/rhamnosus/twntyfr/map_bac
WORKING_DIR=/Volumes/longlunch/seq/LRGC/ruth_meta/data
# Location of the fastq files to map
SAMPLE_DIR=/Volumes/longlunch/seq/LRGC/ruth_meta/gloorg_fastq
# Location of the fasta for the reference index
REFSEQS=/Volumes/longlunch/seq/LRGC/ruth_meta/bowtie2_index/all_genus_orfs_clustered_at_99_unique.fa
OUTPUTFOLDER=/Volumes/longlunch/seq/LRGC/ruth_meta/data/
# Name of your bowtie index. Will output to your working dir
IDX="/Volumes/longlunch/seq/LRGC/ruth_meta/bowtie2_index/nafld"

BIN=/Volumes/longlunch/seq/LRGC/ruth_meta/bin
#----------------------------------------------------------------------------
# Get bowtie version

echo "# Mapper version:"
bowtie2 --version

#/Groups/twntyfr/bin/bowtie2-2.1.0/bowtie2-align version 2.1.0
#64-bit
#Built on ifx5.ebalto.jhmi.edu
#Wed Feb 20 11:34:34 EST 2013
#Compiler: gcc version 4.2.1 (Apple Inc. build 5666) (dot 3)
#Options: -O3 -m64 -msse2 -funroll-loops -g3 
#Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}

# Go to working dir
cd $WORKING_DIR

# Make reference index
if [ -e $IDX.1.bt2 ]; then
	echo -e "\n## Index exists : $IDX\n"
else
	echo -e "\n## Building index : $IDX\n"
	bowtie2-build -f $REFSEQS $IDX > bowtie2-build_log_`date +"%Y-%m-%d_%H%M%S"`.txt
fi


for f in $( ls $SAMPLE_DIR ); do
#	echo "#Working on $f"
	#get the sample name. No idea why I have to split twice
	IFS='\.' read -a array <<< $f
	IFS=' ' read -a array2 <<< $array
	
	SAMPLE=${array2[0]}
	echo "# Working on sample: $SAMPLE"
	DIR="${OUTPUTFOLDER}${SAMPLE}_map"
	
	if [ -d $DIR ]; then
		echo "# Output directory: $DIR exists"
	else
		echo "# Making output directory: $DIR"
		mkdir $DIR
	fi
	
	echo "# Starting mapping on $SAMPLE : `date`"

	bowtie2 -x $IDX -U $SAMPLE_DIR/$f -S $DIR/$SAMPLE.sam -p 10 -N 1 -D 20 -R 3 -L 20 #2> $DIR/errlog.txt

	echo -e "# Done mapping $SAMPLE : `date`\n"

done
    echo -e "\n## All mapping complete\n"

#-------------------------
# Separate loops for the mapping and the partitioning in case one fails
# JM - need to clean these up. Move Perl scripts


for f in $( ls $SAMPLE_DIR ); do
#	echo "#Working on $f"
#get the sample name. No idea why I have to split twice
IFS='\.' read -a array <<< $f
IFS=' ' read -a array2 <<< $array

SAMPLE=${array2[0]}
#	echo "# Working on $SAMPLE"
	DIR="${OUTPUTFOLDER}/${SAMPLE}_map"

    echo "# Making readcounts table for sample: $SAMPLE"
	
	if [ -e "${OUTPUTFOLDER}/headers.txt" ] && [ -s "${OUTPUTFOLDER}/headers.txt" ]; then
		echo "# Headers list exists"
	else
		grep '^@' $DIR/$SAMPLE.sam > "${OUTPUTFOLDER}/headers.txt"
	fi

	if [ -e $DIR/best_hits_${SAMPLE}.out ]; then
		echo "# Read partitioning already complete"
	else
		echo "# Partitioning the reads by genomic position : `date +"%T"`"
#JM		# This prints every single reference (from headers.txt) and counts the fwd and rev reads mapped
#JM		# Probably don't need this. We just need number of reads/ref that are MAPPED
		$BIN/count_reads_SAM.pl $DIR/$SAMPLE.sam > $DIR/best_hits_${SAMPLE}.out
	fi
	
		echo "# Making read count table : `date +"%T"`"

    if [ ! -e $DIR/${SAMPLE}_CDS_counts.txt ]; then
        #Get read count files per gene
        $BIN/count_reads_per_CDS.pl $DIR/best_hits_${SAMPLE}.out > $DIR/${SAMPLE}_CDS_counts_temp.txt
            rm $DIR/best_hits_${SAMPLE}.out
        # Sum up fwd and rev into total reads
        awk -F'\t' 'BEGIN {OFS = "\t"} {s1=$2+$3;print $1,$2,$3,s1}' $DIR/${SAMPLE}_CDS_counts_temp.txt > $DIR/temp.txt
            rm $DIR/${SAMPLE}_CDS_counts_temp.txt
        # Add read length information from headers.txt
        $BIN/add_CDS_length.pl $WORKING_DIR/headers.txt $DIR/temp.txt > $DIR/${SAMPLE}_CDS_counts.txt
            rm $DIR/temp.txt
    else
        echo "# ${SAMPLE}_CDS_counts.txt already made : `date +"%T"`"
    fi

		echo -e "# Read count table for $SAMPLE complete : `date`\n"
done
        echo "## done ALL!"

#---------------------------------------------------------------------------------------------------------
# Notes
#-------------
#bowtie2 flags:

#-p threads
#-un file for unmapped reads (--un $DIR/${SAMPLE}_unmapped.fastq)
#-x the index
#-S name of sam file output
#-U your fastq sample (default format)
#-N 1 number of mismatches (max is 1)
#-D
#-R
#-L

