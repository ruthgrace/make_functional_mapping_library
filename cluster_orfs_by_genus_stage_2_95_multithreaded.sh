for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	echo $directory
	echo $genus
	Rscript sort_seq_by_median_length.r $directory"/"$genus"_multithreaded_cd_hit.txt" $directory"/"$genus"_clustered_at_95_sorted_by_median_length.fa"
	cd-hit-openmp/cd-hit -i $directory"/"$genus"_clustered_at_95_sorted_by_median_length.fa" -o $directory"/"$genus"_95_cd_hit.txt" -c 0.95 -n 5 -T 10
done
