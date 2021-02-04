#!usr/bin/perl
use warnings;

#my $infile = "test";
my $srpssm = $ARGV[0];
my $desp = $ARGV[1];
my $job = $ARGV[2];
my $infile = "$job/$desp.snp.ref.mut.seq";

#my $srpssm = 'SRSF1.pssm';
my %sr = ();

open (IN, $srpssm)||die "$!";
my $head = <IN>;
my $pos = 0;
my $motiflen = 0;
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/);
	$motiflen = $#line;
	for my $i(1..$#line)
	{
		$sr{$i}{$line[0]} = $line[$i];
	}
	$pos +=1;
}
close IN;

my @inname = split(/[\.\/]/,$srpssm);
open (IN, $infile)||die "$!";
open (OUT, ">$infile.$inname[$#inname-1].score")||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\s+/, $_);
	$line[3] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	my $ref_ui = $1;
	my $ref_di = $3;
	my $ref_exon = $2;
	$line[4] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	my $mut_exon = $2;
	my $mut_ui = $1;
	my $mut_di = $3;
	
	my $ref_seq = "";
	my $mut_seq = "";
	
	my @snp = split(/:/, $line[2]);
	my @range = split(/[:-]/, $line[0]);
	my @ref_score = ();
	my @mut_score = ();
	
	for my $sp (0..(length($line[3])- $motiflen))
	{
		my $offs = $sp;
		#if($range[-1] eq '+')
		#{
		#	#$offs = $snp[1] - $range[1] - 300  - $motiflen +1 + $sp;
		#	$offs = $snp[1] - $range[1]   - $motiflen +1 + $sp;
		#}else
		#{
		#	#$offs = $range[2] - 300 - $snp[1]  - $motiflen +1 + $sp;
		#	$offs = $range[2]  - $snp[1]  - $motiflen +1 + $sp;
		#}
		$ref_seq = substr($line[3], $offs , $motiflen); 
		$ref_seq = uc($ref_seq);
		$ref_score[$sp] = PSSM(\$ref_seq, \%sr);
		
		$mut_seq = substr($line[4], $offs , $motiflen); 
		$mut_seq = uc($mut_seq);
		$mut_score[$sp] = PSSM(\$mut_seq, \%sr);		
	}

	my $max_ref_score = (sort {$a <=> $b} @ref_score)[-1];
	my $max_mut_score = (sort {$a <=> $b} @mut_score)[-1];
	print OUT "$line[0]\t$line[1]\t$line[2]\t$max_ref_score\t$max_mut_score\n";
	#print OUT "$line[0]\t$line[1]\t$asscore\t$dsscore\n";
	
}
close IN;
close OUT;

 %sr = ();

sub PSSM
{
	my $tmp  = shift;
	my $pssm = shift;
	
	my @seq = split(//,$$tmp);
	my $jc = 0;
	foreach my $s( sort {$a <=> $b} keys %{$pssm})
	{
		$jc += $$pssm{$s}{$seq[$s-1]};
	}
	
	return $jc;

}
















