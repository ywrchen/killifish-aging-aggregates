# Get the tissue-wise lists of all the genes for GSEA. Values are FC * -log10palue. This is not sorted but will be done later.
# This is modified from previous scripts because now all YC files have the same 4 columns.
# Remember for this script to work directory structure should be ready, which I made manually
use strict;
use warnings;

# ========================================== AGGREGATION PROPENSITY OR AGG OR TL =================================
my $geneIndex = 2;   # Column with gene name 0 based index
my $qvalIndex = 6;   # Column that has p or q value
my $tissueIndex = 0; # Column that has tisse names
my $fcIndex = 5;     # Index of column that has Fold change
my $path = "/Users/yiwenchen/Dropbox/Killifish-Collaboration/Systems_Aggregate_Paper/github/killifish-aging-aggregates/MSAnalysis/SigChanges"; # Path to input files. Change it to your dropbox path to run.
# ==========================================  DATASET =============================================================
# CHOSE THE DATASET AND MATCH THE COLUMN INDEX ABOVE

# I will use these 2 arrays to generate combinations for each comparison and category and separate
# them into different lists in different output directories
my @Comparison = ('OvY', 'TvY', 'TvO');
my @Category = ('AGG', 'Prop', 'TL');

foreach my $set (@Comparison){ # For each comparison
	
	print "> Processing $set\n";
	foreach my $comp (@Category){ # for each category or dataset
		
		my $infile = "$comp\_$set\.csv"; # Read the input file
		print "  $comp $infile\n";
		
		open FILE, "$path/$infile" or die $!;
		my @file = <FILE>;
		shift @file;

		my %SortedGenes; 
		
		# For each protein in the input file do -log10pvalue * FC; sort them for each tissue and put them in SortedGenes.
		foreach (@file){
	
			my @line = split '\,', $_;
			map {$_=~s/\n|\"//g} @line;

			if ($line[0] ne '' && $line[1] ne '' && $line[2] ne '' && $line[3] ne ''){
				$SortedGenes{$line[$tissueIndex]}{$line[$geneIndex]} = -log10($line[$qvalIndex]) * $line[$fcIndex];
			}
			else {
				#print join "\t", @line, "\n";
			}
		}
		#print $SortedGenes{"Gut"}{"acads"};
		
		# For each tissue, generate an output file and print the gene and score
		foreach my $tissue (keys %SortedGenes){

			print "  $tissue\t$set/$comp/$tissue\_$set\_$comp\.csv";
			open (GSEA, ">", "GSEA/GeneSets/$comp/$set/$tissue\_$set\_$comp\.csv") or die $!;
			print GSEA "Protein,mlog10QvalxFC\n";

			print scalar keys %{$SortedGenes{$tissue}}, "\n";

			foreach my $gene (keys %{$SortedGenes{$tissue}}){
		
				print GSEA "$gene,$SortedGenes{$tissue}{$gene}\n";
			}
			close (GSEA);
            #exit; # TEST: This will exit the script and run for only one combination. Uncomment it to run for all combinations.
		}
	}
}

sub log10 {
	my $n = shift;
    return log($n)/log(10);
}

print `date`;
