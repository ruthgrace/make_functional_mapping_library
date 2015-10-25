#!/usr/bin/perl -w
use strict;

#output ORF DNA sequences from a gtf file and a contig fasta file


################
# Need gtf file (download separately), and the "contigs" file (downloaded in package)
#
################

my $contig_file = $ARGV[0] if $ARGV[0];
my $gff_file = $ARGV[1] if $ARGV[1];
print "contig file as first argument gtf file as second argument\n" if !$ARGV[0];
print "contig file as first argument gtf file as second argument\n" if !$ARGV[1];
exit if !$ARGV[0];
exit if !$ARGV[1];

my %contig;
my $lastID;
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $contig_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		$lastID = $l if $l =~ /^>/;
		$contig{$lastID} .= $l if $l !~ /^>/;
	}	
close IN;

#open the gtf file second
my $i = 1;
open (IN, "< $gff_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l !~ /^#/){
			my @l = split/\t/, $l;
			my $id = ">" . $l[0];
			my $start = $l[3] -1;
			my $len = $l[4] - $l[3] -1;
			my $orf = uc( substr($contig{$id}, $start, $len) );
			
			$orf = rev($orf) if $l[6] eq "-";
			
			print "$id|$l[-1]\n$orf\n";
			$i++;
		}
	}	
close IN;


sub rev{
       my $seq;
       my @seq = split//, $_[0];
       for (my $i = @seq - 1; $i >= 0; $i--){
               $seq[$i] =~ tr/ACGT/TGCA/;
               $seq .= $seq[$i];
       }
       return $seq;
}