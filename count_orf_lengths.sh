for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	echo $directory"/"$genus"_all_orfs.fa"
	awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' $directory"/"$genus"_all_orfs.fa"
done
