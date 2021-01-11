# To get other type of ids for the gene names for Yiwen
use strict;
use warnings;

# Outfile 
open OUT, '>Ka_Ks_values_Nfur-to-otherfish_with-XM-XP-Ids.txt' or die $!;

# I use the longest proitein for my analysis. The info is there in the fasta header, so I will read the fasta file and parse the XM and XP ids
local $/ = '>';
open FASTA, 'nfurzeri_LongestProteins_ncbi.txt' or die $!;
my @seqs = <FASTA>;
shift @seqs;

my %SeqInfo;
foreach (@seqs){
	
	# split the l;ines and get the header in the hash
	my @lines = split "\n", $_;
	my $head = shift @lines;
	my @info = split '\|', $head;
	
	$SeqInfo{$info[0]} = "$info[1]\t$info[2]";
}
print "Total ids = ", scalar keys %SeqInfo,"\n";

# Read the file having names and final values- and add the other two ids from SeqInfo hash
local $/ = "\n";
open FH, 'Ka_Ks_values_Nfur-to-otherfish.txt' or die $!;
my @file = <FH>;
my $head = shift @file;
my @head = split "\t", $head;

# Print the head line
print OUT "$head[0]	XM id	XP id	";
print OUT join ("\t", @head[1..$#head]);

# Foreach line, get the id and print in file
foreach (@file){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	if (exists $SeqInfo{$line[0]}){ # If the gene name matches -> print the info and add two ids
		
		print OUT "$line[0]	$SeqInfo{$line[0]}	";
		print OUT join ("\t", @line[1..$#line]), "\n";
	}
	else { # Else -> die
		die "$line[0] id not found. Check!!!";
	}
}

print `date`;

