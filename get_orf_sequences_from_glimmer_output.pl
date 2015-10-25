#!/usr/bin/perl -w
use strict;
use Cwd;
use warnings;

#output ORF DNA sequences from a gtf file and a contig fasta file


################
# Need gtf file (download separately), and the "contigs" file (downloaded in package)
#
################

my $contig_file = $ARGV[0] if $ARGV[0];
my $coord_file = $ARGV[1] if $ARGV[1];
my $outfile = $ARGV[2] if $ARGV[2];
print "contig file as 1st argument, coord file from glimmer as 2nd argument, output file as 3rd argmument \n" if !$ARGV[0];
print "contig file as 1st argument, coord file from glimmer as 2nd argument, output file as 3rd argmument \n" if !$ARGV[1];
print "contig file as 1st argument, coord file from glimmer as 2nd argument, output file as 3rd argmument \n" if !$ARGV[2];
exit if !$ARGV[0];
exit if !$ARGV[1];
exit if !$ARGV[2];

my $dir = getcwd;

unless(open FILE, '>', $outfile) {
    die "\nUnable to open $outfile: $!\n";
}

my %contig;
my $lastID;
my $firstID="";
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $contig_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
			$lastID = $l;
			#for some reason the first ID from the glimmer output is always missing, so save it here
			if ($firstID eq "") {
				$firstID=$l;
			}
		}
		$contig{$lastID} .= $l if $l !~ /^>/;
	}	
close IN;

#open the glimmer coord file second
my $i = 1;
$lastID=$firstID;
open (IN, "< $coord_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
			$lastID = $l;
		}
		if (($l !~ /^>/) && ($l ne "")){
			my @l = split/\s+/, $l;
			my $start = $l[1] -1;
			my $end = $l[2];
			if ($end < $start) {
				my $temp = $start;
				$start = $end;
				$end = $temp;
			}
			my $len = $end - $start +1;
			my $orf = uc( substr($contig{$lastID}, $start, $len) );
			$orf = rev($orf) if $l[3] =~ /^-*/;
			
			print FILE "$lastID|$l[0]\n$orf\n";
			$i++;
		}
	}
close IN;

close FILE;

print "$i ORFS output into $outfile";

sub rev{
       my $seq;
       my @seq = split//, $_[0];
       for (my $i = @seq - 1; $i >= 0; $i--){
               $seq[$i] =~ tr/ACGT/TGCA/;
               $seq .= $seq[$i];
       }
       return $seq;
}