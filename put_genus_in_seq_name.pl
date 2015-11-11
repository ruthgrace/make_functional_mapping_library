#!/usr/bin/perl -w
use strict;
use Cwd;

my $contig_file = $ARGV[0] if $ARGV[0];
my $outfile = $ARGV[1] if $ARGV[1];
my $genus = $ARGV[2] if $ARGV[2];
print "contig file as 1st argument, output file as 2nd argmument, genus as 3rd argument \n" if !$ARGV[0];
print "contig file as 1st argument, output file as 2nd argmument, genus as 3rd argument \n" if !$ARGV[1];
exit if !$ARGV[0];
exit if !$ARGV[1];
exit if !$ARGV[2];

my $dir = getcwd;

unless(open FILE, '>', $outfile) {
    die "\nUnable to open $outfile: $!\n";
}

my %contig;
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $contig_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
      $l =~ s/^>//;
      print FILE ">$genus $l\n";
		}
    else {
      print FILE "$l\n";
    }
	}	
close IN;
close FILE;
