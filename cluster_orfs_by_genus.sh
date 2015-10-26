for directory in `find ./data/orfs -type d`
do
	genus=${directory:12}
	echo $directory
	echo $genus
	usearch -cluster_fast $directory"/"$genus"_all_orfs.fa" -id 1 -centroids $directory"/"$genus"_centroids.ffn" -uc $directory"/"$genus"_clusters_id100.uc"
done

#usearch -cluster_fast data/orfs/Gamma/Gamma_all_orfs.fa -id 1 -centroids ./Gamma_centroids.ffn -uc ./Gamma_clusters_id100.uc
