# To parse paml file and get the dN dS values
use strict;
use warnings;
use List::Util qw(sum);

# open outfile
open OUT, '>Ka_Ks_values_Nfur-to-otherfish.txt' or die $!;
print OUT "Gene	mean Ka/Ks	mean Ka	mean ks	";
print OUT "medaka Ka/Ks	medaka Ka	medaka ks	";
print OUT "Stickleback Ka/Ks	Stickleback Ka	Stickleback ks	";
print OUT "tetraodon Ka/Ks	tetraodon Ka	tetraodon ks	";
print OUT "fugu Ka/Ks	fugu Ka	fugu ks	";
print OUT "zebrafish Ka/Ks	zebrafish Ka	zebrafish ks\n";

# open file to get the gene name
open FH, 'Guidance_joining-order.txt' or die $!;
my %Genes;
foreach (<FH>){
	
	my @line = split "\t", $_;
	map {$_=~s/\n//g} @line;
	
	$Genes{$line[0]} = $line[1];
}

# Open YN00 outfile
local $/ = 'Data set ';
open LIST, "Guidance_YN00_Out.txt" or die $!;
#open OUT, ">dN_dS_values.txt" or die $!;

my @list = <LIST>;
shift @list;
my $i = 1;
#print "*$list[0]*\n";

my $dataset = 1; # To get the gene name from the other file
foreach (@list){ # ----------------------- foreach gene ---------------------------------------------

# Parse the names
my %Names;
while ($_=~/(\d+) \((.+)\) vs. (\d+) \((.+)\)/g){
	
	$Names{$1} = $2;
	$Names{$3} = $4;
	
}

#foreach (keys %Names){
#	print "$_\t$Names{$_}\n";
#}
	
$_ =~/\(B\) Yang & Nielsen \(2000\) method

Yang Z, Nielsen R \(2000\) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17\:32-43

\(equal weighting of pathways\)

(.+)

\(C\) LWL85, LPB93 & LWLm methods
/gs;

my $ynvalues = $1;

my @lines = split "\n", $ynvalues;

my %KaKs;
foreach (@lines){ # ------------- foreach organism --------------
	
	if (($_ !~/^$/g) && ($_ !~/^seq\./g)){

		my @line = split " ", $_;
#		print "$line[0]	$Names{$line[0]}	$Names{$line[1]}	$line[6]	$line[7]	$line[10]\n";

		$KaKs{$Names{$line[0]}}{$Names{$line[1]}} = [$line[6], $line[7], $line[10]];
		$KaKs{$Names{$line[1]}}{$Names{$line[0]}} = [$line[6], $line[7], $line[10]];
	}
} # ----------------------------- end organism ------------------

my @kaks; my @ka; my @ks;
foreach (keys %{$KaKs{'TurKill'}}){
	
#	print "${$KaKs{'TurKill'}{$_}}[0]	${$KaKs{'TurKill'}{$_}}[1]	${$KaKs{'TurKill'}{$_}}[2]\n";
	push @kaks, ${$KaKs{'TurKill'}{$_}}[0];
	push @ka, ${$KaKs{'TurKill'}{$_}}[1];
	push @ks, ${$KaKs{'TurKill'}{$_}}[2];
}

print "$Genes{$dataset}\n";

print OUT "$Genes{$dataset}	";
print OUT mean(@kaks), "\t",mean(@ka), "\t",mean(@ks), "\t" ;
print OUT join ("\t", @{$KaKs{'TurKill'}{'Medaka'}}), "\t";
print OUT join ("\t", @{$KaKs{'TurKill'}{'Sticklebac'}}), "\t";
print OUT join ("\t", @{$KaKs{'TurKill'}{'Tetraodon'}}), "\t";
print OUT join ("\t", @{$KaKs{'TurKill'}{'Fugu'}}), "\t";
print OUT join ("\t", @{$KaKs{'TurKill'}{'Zebrafish'}}), "\n";

$dataset++;

} # -------------------------------------  end gene ------------------------------------------------

print `date`;

sub mean { return @_ ? sum(@_) / @_ : 0 }
