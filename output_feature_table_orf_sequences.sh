for f in ./data/genomes/*.txt
do
	filename=${f:15}
	filename=${filename%_feature_table.txt}
	line=$(head -n 1 "./data/genomes/"$filename"_genomic.fna")
	genus=`echo $line | cut -d " " -f 2`
	echo $f
	echo $filename
	echo $genus
	mkdir -p ./data/orfs/$genus
	./get_orf_sequences_from_feature_table.pl "data/genomes/"$filename"_genomic.fna" "data/genomes/"$filename"_feature_table.txt" "data/orfs/"$genus"/"$filename"_orfs.fa"
done
