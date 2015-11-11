for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	echo $directory
	echo $genus
	./put_genus_in_seq_name.pl $directory"/"$genus"_90_cd_hit.txt" $directory"/"$genus"_90_cd_hit_genus_id.txt" $genus
done
