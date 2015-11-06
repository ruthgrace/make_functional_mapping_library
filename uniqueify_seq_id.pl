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

my $counter = 1;

open (IN, "< $infile") or die "$!\n";
unless(open OUT, '>', $outfile) {
    print "\nUnable to open $outfile: $!\n";
}
	while(defined(my $l = <IN>)){
		chomp $l;
    if ($l =~ /^>/) {
      # strip the >
      $l =~ s/^.//;
      #trim left whitespace
      $l =~ s/^\s+//;
      print OUT "> $counter $l\n";
      $counter++;
    }
    else {
      print OUT "$l\n";
    }
	}	
close IN;
print "Uniqueified $counter sequences\n"
