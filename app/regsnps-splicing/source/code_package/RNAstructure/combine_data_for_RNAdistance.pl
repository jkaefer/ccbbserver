#!usr/bin/perl
use warnings;

#my $infl1 = "HGMD_PRO_INDELS_DM_JAN_12.tsv.raw.cor.bed.seq.rnafold"; # ref sequence 2nd structure
#my $infl2 = "HGMD_PRO_INDELS_DM_JAN_12.tsv.mut.cor.bed.seq.use.rnafold"; # mut sequence 2nd structure
my $infl1 = $ARGV[0];
my $infl2 = $ARGV[1];

my @ref = ();
my @mut = ();

my @seq = ();
my $count = 0;
open (IN, "$infl1") || die "$!";
while(<IN>)
{
	$count++;
	if ($count < 3)
	{
		push @seq, $_;
	}
	else
	{
		push @seq, $_;
		push @ref, [@seq];
		@seq = ();
		$count = 0;
	}
}
close IN;


open (IN,"$infl2")||die "$!";
while(<IN>)
{
	$count++;
	if ($count < 3)
	{
		push @seq, $_;
	}
	else
	{
		push @seq, $_;
		push @mut, [@seq];
		@seq = ();
		$count = 0;
	}
}
close IN;

for my $i (0..$#ref)
{
	foreach my $j ( @{$ref[$i]} )
	{
		print "$j";
	}
	foreach my $k ( @{$mut[$i]} )
	{
		print "$k";
	}
}
