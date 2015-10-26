for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	echo $directory
	echo $genus
	cat $directory/* > $directory"/"$genus"_all_orfs.fa"
done
