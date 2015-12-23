from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys

#this script takes two arguments: the input file followed by the output file.
input_file = sys.argv[1]
output_file = sys.argv[2]

handle = open(input_file, "rU")
sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

keys = sequences.keys()

output_handle = open(output_file, "w")
for key in keys:
  coding_dna = sequences[key].seq
  sequences[key].seq = coding_dna.translate(table=11)
  SeqIO.write(sequences[key], output_handle, "fasta")
  print("processed seq with key " + key)

output_handle.close()
print("done")
