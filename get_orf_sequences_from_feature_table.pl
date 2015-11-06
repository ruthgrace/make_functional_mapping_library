#!/usr/bin/perl -w
use strict;
use Cwd;

my $contig_file = $ARGV[0] if $ARGV[0];
my $feature_table = $ARGV[1] if $ARGV[1];
my $outfile = $ARGV[2] if $ARGV[2];
print "contig file as 1st argument, feature table file as 2nd argument, output file as 3rd argmument \n" if !$ARGV[0];
print "contig file as 1st argument, feature table file as 2nd argument, output file as 3rd argmument \n" if !$ARGV[1];
print "contig file as 1st argument, feature table file as 2nd argument, output file as 3rd argmument \n" if !$ARGV[2];
exit if !$ARGV[0];
exit if !$ARGV[1];
exit if !$ARGV[2];

my $dir = getcwd;

unless(open FILE, '>', $outfile) {
    die "\nUnable to open $outfile: $!\n";
}

my %contig;
my $lastID;
#open the contig file first and save string in a hash keyed by contig name
open (IN, "< $contig_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l =~ /^>/) {
			$lastID = $l;
			($lastID) = $lastID =~ /\A([^:\s]+)/;
		}
		$contig{$lastID} .= $l if $l !~ /^>/;
	}	
close IN;

#open the feature table file second
my $i = 1;
open (IN, "< $feature_table") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
		if ($l !~ /^#/){
			my @l = split/\t/, $l;
			if ($l[0] eq "CDS") {
				my $id = ">" . $l[6];
				my $start = $l[7] -1;
				my $len = $l[8] - $l[7] +1;
				my $orf = uc( substr($contig{$id}, $start, $len) );
				
				$orf = rev($orf) if $l[9] eq "-";
				
				print FILE "$id|$l[10]\n$orf\n";
				$i++;
			}
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