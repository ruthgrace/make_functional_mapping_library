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

my %seqs;
my %lengths;
my $lastID = "";
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $infile") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
      if ($lastID != "") {
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
foreach my $id (sort { $lengths{$a} <=> $lengths{$b} } keys %lengths) {
  print "Printing sequence $id with length $lengths{$id}\n";
  print OUT "$id\n";
  print OUT "$seqs{$id}";
}
close OUT;
