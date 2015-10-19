fna='.fna'
prefix='iterated_glimmer_annotation_'
for f in ./data/genomes/*.gff
do
	glimmer/glimmer3.02/scripts/g3-iterated.csh ${f%.gff}$fna $prefix${f%.gff}
done
