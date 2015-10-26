for f in ./data/genomes/*.gff
do
	filename=${f:15}
	filename=${filename%_genomic.gff}
	line=$(head -n 1 "./data/genomes/"$filename"_genomic.fna")
	genus=`echo $line | cut -d " " -f 2`
	if  [[ $genus == \[* ]] ;
	then
	    len=${#genus}
	    genus=${genus:1:(len-2)}
	fi
	echo $f
	echo $filename
	echo $genus
	mkdir -p ./data/orfs/$genus
	./get_orf_sequences_from_glimmer_output.pl "data/genomes/"$filename"_genomic.fna" "data/glimmer_output/"$filename"_genomic.coords" "data/orfs/"$genus"/"$filename"_orfs.fa"
done
