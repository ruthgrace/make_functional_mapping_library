fna='.fna'
prefix='iterated_glimmer_annotation_'
for f in ./data/genomes/*.gff
do
  filename=${f:15}
  filename=${filename%.gff}
  output=./data/glimmer_output/$filename
  seq=${f%.gff}
  seq=$seq.fna
  echo $f
  echo $filename
  echo $output
  echo $seq
  glimmer/glimmer3.02/bin/long-orfs -n -t 1.15 $seq $output.longorfs
  glimmer/glimmer3.02/bin/extract -t $seq $output.longorfs > $output.train
  glimmer/glimmer3.02/bin/build-icm -r $output.icm < $output.train
  glimmer/glimmer3.02/bin/glimmer3 -o50 -g110 -t30 $seq $output.icm $output.run1
  tail +2 $output.run1.predict > $output.coords
done
