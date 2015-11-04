for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	awk -v path=$directory"/"$genus"_all_orfs.fa" 'BEGIN {ORS="\t"} /^>/ {if (seqlen){printf "%s\n", seqlen}; print path;print;seqlen=0;next; } { seqlen = seqlen +length($0)}END{printf "%s\n", seqlen}' $directory"/"$genus"_all_orfs.fa"
done
