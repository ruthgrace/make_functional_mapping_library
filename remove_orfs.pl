#!/usr/bin/perl -w
use strict;
use Cwd;
use warnings;

#output ORF DNA sequences from a gtf file and a contig fasta file


################
# Need gtf file (download separately), and the "contigs" file (downloaded in package)
#
################

my $unwanted_orf_file = $ARGV[0] if $ARGV[0];
print "This function takes in one argument: the path of a table with columns filepath, sequence name, and sequence length, for each sequence that should be removed\n" if !$ARGV[0];
exit if !$ARGV[0];

my %unwanted_seqs;

open (IN, "< $unwanted_orf_file") or die "$!\n";
	while(defined(my $l = <IN>)){
		chomp $l;
    my @items = split /\s/, $l;
    my $filename = $items[0];
    my $seqlength = $items[-1];
    # get sequence identifier
    $l =~ s/(^\Q$filename\E)|(\Q$seqlength\E$)//g;
		# trim sequence identifier
		$l =~ s/^\s+|\s+$//g;
		print "saving filename $filename sequence id $l into unwanted_seqs\n";
    $unwanted_seqs{$filename}{$l} = 1;
	}	
close IN;

my $keep_seq = 1;

my $orf_path = "./data/orfs/";
my $orf_file_ending = "_all_orfs.fa";
my $filtered_orf_file_ending = "_filtered_orfs.fa";

my @dirs = grep { -d } glob "$orf_path*";
my $genus_folder;
my $genus;
foreach $genus_folder (@dirs) {
  $genus = (split(/\//, $genus_folder))[-1];
  my $filename = "$genus_folder/$genus$orf_file_ending";
  my $outfilename = "$genus_folder/$genus$filtered_orf_file_ending";
  unless(open FILE, '>', $outfilename) {
      print "\nUnable to open $outfilename: $!\n";
  }
  else {
    unless(open IN, '<', $filename) {
        print "\nUnable to open $filename: $!\n";
    }
    else {
    	while(defined(my $l = <IN>)){
				# trim sequence identifier
    		$l =~ s/^\s+|\s+$//g;
    		if ($l =~ /^>/) {
					print "Checking filename $filename sequence id $l for existence in unwanted seqs\n";
          if (exists $unwanted_seqs{$filename}{$l}) {
						print "Sequence found\n";
            $keep_seq = 0;
          }
          else {
            $keep_seq = 1;
            print FILE "$l\n";
          }
    		}
        elsif ($keep_seq==1) {
          print FILE "$l\n";
        }
    	}
      close IN;
      print "Filtered orfs written to $outfilename\n";
    }
  }
}
