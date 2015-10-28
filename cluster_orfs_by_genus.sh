for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	echo $directory
	echo $genus
	Rscript sort_seq_by_median_length.r $directory"/"$genus"_all_orfs.fa" $directory"/"$genus"_all_orfs_sorted_by_median_length.fa"
	cd-hit-v4.6.1-2012-08-27/cd-hit -i $directory"/"$genus"_all_orfs_sorted_by_median_length.fa" -o $directory"/"$genus"_cd_hit.txt" -c 1.00 -n 5
done
