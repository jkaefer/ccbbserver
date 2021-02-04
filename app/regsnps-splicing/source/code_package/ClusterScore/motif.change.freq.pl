#usr/bin/perl
use warnings;

## calculate the normalized motif/hexamer gain/loss freq, and affected location freq
my $pref = $ARGV[0];#'hgmd.exonic' or 'neutral'
my $infile = "/N/slate/jkaefer/splicingdb/motif/known.predicted.motif";
my %known_motif = ();
my %hexamer = ();
my %motif_count = ();
my %motif_affected_loc = ();
my $pseudo_count = 1;
my $total_count = 0; # for motif gain/loss, not location
open(IN, $infile)||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	my @tmp = split(/\./regsnps-splicing/, $line[0]);
	my $id = join(".", @tmp[0..($#tmp-1)]);
	$line[1] = uc($line[1]);
	if ($line[0] =~ /chasin|burge/){
		my @hmer = split(/\|/, $line[1]);
		foreach my $h (@hmer){
			$motif_count{"$h"}{'gain'} =  $pseudo_count;
            $motif_count{"$h"}{'loss'} =  $pseudo_count;
            $total_count +=2;
			for my $i (1..length($h)){
				$motif_affected_loc{$h}{'gain'}{$i} =1;
				$motif_affected_loc{$h}{'loss'}{$i} =1;
			}
		}
	}else{
		$motif_count{$id}{'gain'} = $pseudo_count;
        $motif_count{$id}{'loss'} = $pseudo_count;
        $total_count +=2;
		my @group = split(/\|/, $line[1]);
		for my $i (1..length($group[0])){
			$motif_affected_loc{$id}{'gain'}{$i} =1;
			$motif_affected_loc{$id}{'loss'}{$i} =1;
		}
	}
	my @hex = split(/\|/, $line[1]);
	foreach my $h (@hex){
		$hexamer{$h} = $id;
	}
}
close IN;
$total_count -= 98;

open (IN, "$pref.snp.ref.mut.seq.ESR.gain-loss.count")||die "$!";
while(<IN>)
{
    s/\s+$//;
    my @line = split(/\t/, $_);
    my @info = split(/:/, $line[0]);
    my $typ  = '';## ese:gain, ese:loss, ess:gain, ess:loss
    $motif_count{$info[1]}{$info[0]} += $line[1];
    $total_count += $line[1];
    #print "$motif_count{$info[1]}{$info[0]}\n";
}
close IN;

open (OUT, ">$pref.gain-loss.motif.freq")||die "$!";
foreach my $m (keys %motif_count)
{
	#if ($motif_count{$m}{'gain'} == 1){
	#	$motif_count{$m}{'gain'} = 0;
	#}
    $motif_count{$m}{'gain'} /= $total_count;
    #print "$m";
	#if ($motif_count{$m}{'loss'} ==1 ){
	#	$motif_count{$m}{'loss'} =0;
	#}
    $motif_count{$m}{'loss'} /= $total_count;
    #print "$motif_count{$m}{'loss'}\n";
	#my $ratio = sprintf("%.3f", $motif_count{$m}{'gain'}/$motif_count{$m}{'loss'});
	print OUT "$m\t$motif_count{$m}{'gain'}\t$motif_count{$m}{'loss'}\n";
}
close OUT;

open (GA, "$pref.snp.ref.mut.seq.ESR.gain.location.count")||die "$!";
while(<GA>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	$motif_affected_loc{$line[0]}{'gain'}{$line[2]} += $line[3];
}
close GA;
open (GA, "$pref.snp.ref.mut.seq.ESR.loss.location.count")||die "$!";
while(<GA>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	$motif_affected_loc{$line[0]}{'loss'}{$line[2]} += $line[3];
}
close GA;

open(OUT, ">$pref.motif.location.affected.freq")||die "$!";
foreach my $m (keys %motif_affected_loc){
	my $tc = 0;
	print OUT "$m";
	foreach my $i (keys %{$motif_affected_loc{$m}{'gain'}}){
		$tc += $motif_affected_loc{$m}{'gain'}{$i};
	}
	foreach my $i (keys %{$motif_affected_loc{$m}{'gain'}}){
		#if ($motif_affected_loc{$m}{'gain'}{$i} ==1){
		#	$motif_affected_loc{$m}{'gain'}{$i} = 0;
		#}
		$motif_affected_loc{$m}{'gain'}{$i} /= $tc;
		my $f = sprintf("%.10f", $motif_affected_loc{$m}{'gain'}{$i});
		print OUT "\t$f";
	}
	
	$tc = 0;
	foreach my $i (keys %{$motif_affected_loc{$m}{'loss'}}){
		$tc += $motif_affected_loc{$m}{'loss'}{$i};
	}
	foreach my $i (keys %{$motif_affected_loc{$m}{'loss'}}){
		#if ($motif_affected_loc{$m}{'loss'}{$i} ==1){
		#	$motif_affected_loc{$m}{'loss'}{$i} = 0;
		#}
		$motif_affected_loc{$m}{'loss'}{$i} /= $tc;
		my $f = sprintf("%.10f",$motif_affected_loc{$m}{'loss'}{$i});
		print OUT "\t$f";
	}
	print OUT "\n";
}
close OUT;
