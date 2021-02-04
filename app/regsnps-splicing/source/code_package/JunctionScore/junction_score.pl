#!usr/bin/perl
use warnings;

my $desp = $ARGV[0];
my $job = $ARGV[1];
my $infile = "$job/$desp.snp.ref.mut.seq";
my $accsite   = "/N/slate/jkaefer/splicingdb/motif/acceptorsite.pssm";
my $donrosite = "/N/slate/jkaefer/splicingdb/motif/donorsite.pssm";

my %accpssm = ();
my %donorpssm = ();
open (IN, $accsite)||die "$!";
my $head = <IN>;
my $pos = 0;
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/);
	$accpssm{$pos}{'G'} = $line[1];
	$accpssm{$pos}{'A'} = $line[2];
	$accpssm{$pos}{'C'} = $line[3];
	$accpssm{$pos}{'T'} = $line[4];
	$pos +=1;
}
close IN;

$pos = 0;
open (IN, $donrosite)||die "$!";
$head = <IN>;
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	$donorpssm{$pos}{'G'} = $line[1];
	$donorpssm{$pos}{'A'} = $line[2];
	$donorpssm{$pos}{'C'} = $line[3];
	$donorpssm{$pos}{'T'} = $line[4];
	$pos +=1;
}
close IN;

open (OUT, ">$infile.junction-score")||die "$!";
open (MUT, ">$infile.junction-score.mut")||die "$!";
open (IN, $infile)||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
#### ref sequence
	$line[3] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	my $exon = $2;
	my $upintron = $1;
	my $dnintron = $3;

	my $jcsite = substr($upintron, length($upintron)-13, 13); ## acceptor site
	$jcsite   .= substr($exon,0,1);
	$jcsite = uc($jcsite);
	my $asscore = junction(\$jcsite, \%accpssm);

	$jcsite  = substr($exon, length($exon)-3,3);#donor site
	$jcsite .= substr($dnintron, 0, 7);
	$jcsite = uc($jcsite);
	my $dsscore = junction(\$jcsite, \%donorpssm);

	print OUT "$line[0]\t$line[1]\t$line[2]\t$asscore\t$dsscore\n";
	
#### mut sequence
	$line[4] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	$exon = $2;
	$upintron = $1;
	$dnintron = $3;

	$jcsite = substr($upintron, length($upintron)-13, 13); ## acceptor site
	$jcsite   .= substr($exon,0,1);
	$jcsite = uc($jcsite);
	$asscore = junction(\$jcsite, \%accpssm);

	$jcsite  = substr($exon, length($exon)-3,3);#donor site
	$jcsite .= substr($dnintron, 0, 7);
	$jcsite = uc($jcsite);
	$dsscore = junction(\$jcsite, \%donorpssm);

	print MUT "$line[0]\t$line[1]\t$line[2]\t$asscore\t$dsscore\n";
	

}
close IN;
close OUT;
close MUT;

%accpssm = ();
%donorpssm = ();


sub junction
{
	my $tmp  = shift;
	my $pssm = shift;
	
	my @seq = split(//,$$tmp);
	my $jc = 0;
	foreach my $s( sort {$a <=> $b} keys %{$pssm})
	{
		$jc += $$pssm{$s}{$seq[$s]};
	}
	
	return $jc;

}





