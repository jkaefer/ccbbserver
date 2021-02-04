#!usr/bin/perl
use warnings;

my $infile = $ARGV[0];
my $refile = $ARGV[1];
open (IN,"$infile")||die "$!";
my $max = 0;
my $tpr = 0;
my $fpr = 0;
my $tnr = 0;
my $fnr = 0;

my $acc = 0;
my $mcc = 0;

my $corr = 0;
my $errc = 0;
my $threshold = -1;
while(<IN>)
{
	s/\s+$//;
	if (!/@/ and $_ ne '')
	{
		my @line = split(/,/, $_);
		my $tp = $line[0];
		my $fn = $line[1];
		my $fp = $line[2];
		my $tn = $line[3];
		
		my $p = $tp + $fn +1;
		my $n = $tn + $fp +1;
		my $pt = $tp + $fp+1;
		my $nt = $tn + $fn+1;


		if ( $p*$n*$pt*$nt != 0)
		{
			$mcc = ($tp*$tn - $fp*$fn)/sqrt($p*$n*$pt*$nt);
			if ($mcc > $max)
			{
				$max = $mcc;
				$fpr = $line[4];
				$tpr = $line[5];
				$tnr = $tn/($tn+$fp);
				$fnr = $fn/($fn+$tp);

				$acc = ($tp+$tn)/($tp+$tn+$fp+$fn);

				$corr = $tp + $tn;
				$errc = $fn + $fp;
				$threshold = $line[$#line];
			}
		}
	}
}
close IN;

my @auc = ();
open (IN, $refile)||die "$!";
while(<IN>)
{
	s/\s+$//;
	if (/^Weighted Avg/)
	{
		my @line = split(/\s+/, $_);
		push @auc, $line[$#line];
	}
}
close IN;

my $roc = 0;
if ($#auc > 0)
{
	$roc = $auc[1];
}else{
	$roc = $auc[0];
}

print "mcc:$max\tfpr:$fpr\ttpr:$tpr\tfnr:$fnr\ttnr:$tnr\tacc:$acc\tcorrect err:$corr,$errc\tAUC:$roc\tthreshold:$threshold\n";
