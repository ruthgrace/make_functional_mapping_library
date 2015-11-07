#!/usr/bin/perl -w
use strict;
use Cwd;
use warnings;

my $infile = $ARGV[0] if $ARGV[0];
my $outfile = $ARGV[1] if $ARGV[1];
print "This function takes in two arguments: the path of the input fasta file, and the path of the output fasta file\n" if !$ARGV[0];
print "This function takes in two arguments: the path of the input fasta file, and the path of the output fasta file\n" if !$ARGV[1];
exit if !$ARGV[0];
exit if !$ARGV[1];

print "Reading sequences\n";

my %seqs;
my %lengths;
my $lastID = "";
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $infile") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
      if ($lastID ne "") {
        $lengths{$lastID} = length $seqs{$lastID};
      }
			$lastID = $l;
      $seqs{$lastID} = "";
		}
		else {
      $seqs{$lastID} .= "$l\n";
    }
	}	
close IN;

unless(open OUT, '>', $outfile) {
    print "\nUnable to open $outfile: $!\n";
}
print "Finding median\n";
# find median
my @seqlengths = sort { $lengths{$a} <=> $lengths{$b} } keys %lengths;
my $num_seqs = int(scalar @seqlengths);
my $median = $lengths{$seqlengths[int( $num_seqs/2 )]};
print "Median length is $median, total number of sequences is $num_seqs.\n";
print "Calculating the difference between each sequence length and the median\n";
# subtract each seqlength by median (absolute)
foreach my $id (keys %lengths) {
  $lengths{$id} = abs($lengths{$id} - $median);
}

#sort and output
foreach my $id (sort { $lengths{$a} <=> $lengths{$b} } keys %lengths) {
  print OUT "$id\n";
  print OUT "$seqs{$id}";
}
close OUT;

print "Done\n";
