#!usr/bin/perl
use warnings;
use List::Util qw[min max];

#my $infile = "hgmd.exonic.snp.ref.mut.seq";
my $infile = $ARGV[0];
my $job = $ARGV[1];

my $outfile = "$job/$infile.efscore.input";
open (OUT, ">$outfile")||die "$!"; 
open(REF, ">$job/$infile.ref.ranfold.input")||die "$!";
open(MUT, ">$job/$infile.mut.ranfold.input")||die "$!";

open (IN, "$job/$infile")||die "$!";

my $lc = 0;

while(<IN>)
{
	$lc +=1;
#	if($lc < 5000){
	s/\s+$//;
	my @line = split(/\t/, $_);
#	$line[3] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
#	$line[4] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	
	my $mask = $line[3] ^ $line[4];
	$mask =~ /[^\0]/g;
	my $snp_pos = $-[0];
	my $affected_region_ref = uc(substr($line[3], $snp_pos-20, 41));
	my $affected_region_mut = uc(substr($line[4], $snp_pos-20, 41));
	print OUT ">$line[0]\t$line[1]\t$line[2]\n$affected_region_ref\n>$line[0]\t$line[1]\t$line[2]\n$affected_region_mut\n";

	$affected_region_ref = uc(substr($line[3], $snp_pos-30, 61));
	$affected_region_mut = uc(substr($line[4], $snp_pos-30, 61));
	print REF ">$line[0]\t$line[1]\t$line[2]\t\n$affected_region_ref\n";
	print MUT ">$line[0]\t$line[1]\t$line[2]\t\n$affected_region_mut\n";
#}
}
close IN;
close OUT;
close REF;
close MUT;



my $efscorefl = "$job/$infile.efscore";
my $status = system("perl  ./regsnps-splicing/source/code_package/EFScore/GetSecondaryStructureValues.perl -f $outfile -l 6 -o $efscorefl -method PU -j $job");

open (OUT, ">$efscorefl.value")||die "$!";
#print OUT "$status\n";
	open(IN, $efscorefl)||die "$!";
	while(<IN>)
	{
		s/\s+$//;
		my @line = split(/\t/, $_);
		my @ref_ef_score= split(/;/, $line[$#line]);
		my $off = 1; ## (51-15)/2 +1
		my $efv = 0;
		#for my $i (1..$#ref_ef_score){
		for my $i (15..21){	
			$efv += $ref_ef_score[$i];
		}
		#$efv /= ($#ref_ef_score);
		$efv /= 7;
		print OUT "$line[0]\t$line[1]\t$efv\t";

		my $al = <IN>;
		$al =~ s/\s+$//;
		@line = split(/\t/, $al);
		my @mut_ef_score = split(/;/, $line[$#line]);
		$efv = 0;
		for my $i (15..21){
			$efv += $mut_ef_score[$i];
		}
		$efv /= 7;

		print OUT "$efv\n";
	}
	close IN;
	close OUT;




