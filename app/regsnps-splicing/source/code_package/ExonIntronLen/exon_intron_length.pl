#!usr/bin/perl
#use warnings;

## get the length of exons which has SNP and flanking introns
my $desp = $ARGV[0];
my $job = $ARGV[1];

my $infile = "$job/$desp.snp.ref.mut.seq";
my %exons = ();
open (IN, $infile)||die "$!";
my $count = 0;
while(<IN>){
	s/\s+$//;
	my @line = split(/\t/, $_);
	$exons{$line[0]} = $count;
	$count ++;
}
close IN;

my %length = ();

open (IN, '/N/slate/jkaefer/splicingdb/hg19.bed')||die "$!";
my $head = <IN>;
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	my @leftb = split(/,/, $line[9]);
	my @rightb = split(/,/, $line[10]);## no matter on forward strand or backward strand, leftb is less than rightb

	my @leftb_w_utr = @leftb;	
	my @rightb_w_utr = @rightb;
	
	my $ssl = 0;
	my $ssr = 0;
# remove UTR region
	for my $i (0..$#leftb)
	{
		if ($leftb[$i] <= $line[6] and $rightb[$i] >= $line[6])
		{
			$ssl = $i;
		}		
		if ($leftb[$i] <= $line[7] and $rightb[$i] >= $line[7])
		{
			$ssr = $i;
		}
	}
	my $loffset = $ssl;
	my $roffset = $#rightb - $ssr;

	@leftb  = @leftb[ $ssl..$ssr];
	@rightb = @rightb[$ssl..$ssr];
	$leftb[0]   = $line[6];
	$rightb[-1] = $line[7];
	
	if ($#leftb == $#rightb)
	{
		for my $i (0..$#leftb)
		{
			my $pl = $leftb[$i] -300+1;
			my $pr = $rightb[$i]+300;
			my $id = "$line[2]:$pl-$pr:$line[3]";
			if ($exons{$id} ne '')
			{
				my $lintronl = -1;
				if($i-1 >=0| $ssl > 0)
				{	
					$lintronl = $leftb[$i] - $rightb[$i-1] +1;
					if($ssl > 0)
					{
						$lintronl = $leftb_w_utr[$i+$loffset] - $rightb_w_utr[$i-1+$loffset] +1;
					}
				}
				my $rintronl = -1;
				if ($i+1 <= $#leftb | $ssr < $#rightb_w_utr)
				{
					$rintronl = $leftb[$i+1] - $rightb[$i] +1;
					if($ssr < $#rightb_w_utr)
					{
						$rintronl = $leftb_w_utr[$i+1+$loffset] - $rightb_w_utr[$i + $loffset ] +1;
					}
				}
				my $exonl = $rightb[$i] - $leftb[$i] +1;
				if ($line[3] eq '-')
				{
					my $tmp = $lintronl;
					$lintronl = $rintronl;
					$rintronl = $tmp;
				}
				$length{$id} = [($lintronl-1,$exonl-1,$rintronl-1)];
			}

		}
		 
	}
}
close IN;

open (IN, $infile)||die "$!";
open (OUT , ">$infile.exon-intron.len")||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	print OUT "$line[0]\t$line[1]\t$line[2]";
	
	if ($length{$line[0]} ne '')
	{
		foreach my $e (@{$length{$line[0]}})
		{
			print OUT "\t$e";
		}
		print OUT "\n";
	}
	else{
		print OUT "\tnone\tnone\tnone\n";
	}
}
close IN;
close OUT;
%exons = ();
