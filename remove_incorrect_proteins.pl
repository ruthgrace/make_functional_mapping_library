#!/usr/bin/perl -w
use strict;
use Cwd;
use warnings;

my $infile = $ARGV[0] if $ARGV[0];
my $outfile = $ARGV[1] if $ARGV[1];
my $badfile = $ARGV[2] if $ARGV[2];
print "This function takes in two arguments: the path of the input fasta file, and the path of the output fasta file, and the path of the incorrectly annotated sequences file (these sequences have - or * in them)\n" if !$ARGV[0];
print "This function takes in two arguments: the path of the input fasta file, and the path of the output fasta file, and the path of the incorrectly annotated sequences file (these sequences have - or * in them)\n" if !$ARGV[1];
print "This function takes in two arguments: the path of the input fasta file, and the path of the output fasta file, and the path of the incorrectly annotated sequences file (these sequences have - or * in them)\n" if !$ARGV[2];
exit if !$ARGV[0];
exit if !$ARGV[1];
exit if !$ARGV[2];

my $asterisk = "*";
my $dash = "-";

my %seqs;
my $lastID = "";
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $infile") or die "$!\n";
unless(open OUT, '>', $outfile) {
    print "\nUnable to open $outfile: $!\n";
}
unless(open BAD, '>', $badfile) {
    print "\nUnable to open $badfile: $!\n";
}
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
			if ($lastID ne "") {
        if ($seqs{$lastID} =~ /.*[*]$/) {
          $seqs{$lastID} =~ s/[*]$//;
        }
        if (index($seqs{$lastID}, $asterisk) != -1 || index($seqs{$lastID}, $dash) != -1) {
  				print BAD "$lastID\n";
  				print BAD "$seqs{$lastID}";
        }
        else {
  				print OUT "$lastID\n";
  				print OUT "$seqs{$lastID}";
        }
			}
			$lastID = $l;
      $seqs{$lastID} = "";
		}
		else {
      $seqs{$lastID} .= "$l\n";
    }
	}	
close IN;
close OUT;
close BAD;
