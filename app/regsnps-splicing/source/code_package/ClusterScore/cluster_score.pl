#!usr/bin/perl
#use warnings;
#calculate the occurence/frequency of each motif in each region ## only in the host exon
## and calcualte the ESR-gain-loss and ESR score

#my $infile = "/data/pro/as/data/motif/known.predicted.motif";
my $infile = "/N/slate/jkaefer/splicingdb/motif/known.predicted.motif";# this file put chasin/burge predicted motif in the first two rows
#my $infile = "test.motif";
#my $regfl = "/data/pro/as/data/hg19.exons.7regions.out";
#my $regfl = "hgmd.exonic.snp.ref.mut.seq";
my $regfl = $ARGV[0];
my $job = $ARGV[1];

my @region = ('AI1', 'A', 'AI2');
my %known_motif = ();
my %sffreq = ();
my %hexamer = ();

open (IN, $infile)||die "$!";
open (OUT, ">$job/$regfl.motif_cluster.score")||die "$!";

while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	my @tmp = split(/\./, $line[0]);
	my $id = join(".", @tmp[0..($#tmp-1)]);
	$line[1] = uc($line[1]);
	if ($line[0] =~ /chasin|burge/){
		my @hmer = split(/\|/, $line[1]);
		foreach my $h (@hmer){
			$known_motif{"$id.$h"} =  $h;
		}
	}else{
		$known_motif{$id} = $line[1];	
	}
	
	my @hex = split(/\|/, $line[1]);
	foreach my $h (@hex){
		$hexamer{$h} = $id;
	}
}
close IN;

print OUT "snp\texonID\tgene";
	#foreach my $pos (@region){
		foreach my $kk (sort keys %known_motif){
			print OUT "\t$kk";
		}
	#}
	#foreach my $pos (@region){
		foreach my $kk (sort keys %known_motif){
			print OUT "\t$kk.mut";
		}
	#}
	print OUT "\n";


open (IN, "$job/$regfl")||die "$!";
#my $head = <IN>;
my %esr_gain_loss = (); # with position
my %gain_loss = ();# without position
my %snp_affected_motif = ();

while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	
	#get the snp position
	hexamer_score(\$line[3], \$line[4],\%hexamer, \%esr_gain_loss, \%gain_loss, \%snp_affected_motif, $line[2]);	
	
	factor_occurence(\@line, \%known_motif , \%sffreq); ## ref sequence
	print OUT "$line[0]\t$line[1]\t$line[2]";
	
	#foreach my $pos (@region){
		foreach my $kk (sort keys %{$sffreq{$line[0]}}){
			print OUT "\t$sffreq{$line[0]}{$kk}{'A'}";
		}		
	#}
	%sffreq = ();
	my @mut_seq = @line[0,1,2,4];
	factor_occurence(\@mut_seq, \%known_motif , \%sffreq);## mutated sequence
	#foreach my $pos (@region){
		foreach my $kk (sort keys %{$sffreq{$line[0]}}){
			print OUT "\t$sffreq{$line[0]}{$kk}{'A'}";
		}		
	#}
		
	print OUT "\n";
	%sffreq = ();
}

close IN;
close OUT;

##### motif-location affected frequency
open (GA, ">$job/$regfl.ESR.gain.location.count")||die "$!";
foreach my $m (sort keys %{$esr_gain_loss{'gain'}}){
	foreach my $l (sort keys %{$esr_gain_loss{'gain'}{$m}}){
		if($m =~ /^[ACGTN]+/){
			print GA "$m\t$hexamer{$m}\t$l\t$esr_gain_loss{'gain'}{$m}{$l}\n";
		}else{
			my @mo = grep {$hexamer{$_} eq $m} keys %hexamer;
			my $m_cluster = join('|', @mo);
			print GA "$m\t$m_cluster\t$l\t$esr_gain_loss{'gain'}{$m}{$l}\n";
		}
		
	}	
}
close GA;

open (LO, ">$job/$regfl.ESR.loss.location.count")||die "$!";
foreach my $m (sort keys %{$esr_gain_loss{'loss'}}){
	foreach my $l (sort keys %{$esr_gain_loss{'loss'}{$m}}){
		if($m =~ /^[ACGTN]+/){
			print LO "$m\t$hexamer{$m}\t$l\t$esr_gain_loss{'loss'}{$m}{$l}\n";
		}else{
			my @mo = grep {$hexamer{$_} eq $m} keys %hexamer;
			my $m_cluster = join('|', @mo);
			print LO "$m\t$m_cluster\t$l\t$esr_gain_loss{'loss'}{$m}{$l}\n";
		}
		
	}
}
close LO;

## motif affected frequency
open (GA, ">$job/$regfl.ESR.gain-loss.count")||die "$!";
foreach my $m (sort keys %gain_loss){
	#foreach my $s (sort keys %{$gain_loss{$m}}){
		#print GA "$m\t$s\t$gain_loss{$m}{$s}\n";
	print GA "$m\t$gain_loss{$m}\n";
	#}
	
}
close GA;

open (GA, ">$job/$regfl.snp-affected-motif")||die "$!";
foreach my $s (sort keys %snp_affected_motif){
     print GA "$s\t$snp_affected_motif{$s}{'gain'}\t$snp_affected_motif{$s}{'loss'}\n";
}
close GA;


%known_motif = ();
%sffreq = ();
%hexamer = ();
%esr_gain_loss = ();
%gain_loss = ();
%snp_affected_motif = ();




sub factor_occurence
{#score calculation in exon	
	my $seq = shift;
	my $mot = shift;
	my $count = shift;
	my @line = @$seq;
	$line[3] =~ /^([acgtn]+)([ACGTN]+)([acgtn]+)$/;
	my $upintron = $1;
	my $exon     = $2;
	my $dnintron = $3;
	$exon = substr($upintron, -30).$exon.substr($dnintron, 0, 30);
	$exon = uc($exon);
	
	foreach my $m (keys %$mot)
	{
		my $count_u = 0;
		my $count_e = 0;
		my $count_d = 0;
		my $total_match = 0;
	
		my @group = split(/\|/, $$mot{$m});
		#@group = ($group[0]);
		my @matches = (0) x length($line[3]);
		my @pos = ();
		my $gap = 0;
		foreach  my $g (@group){
			#$g = "T.TTT";
			$gap = length($g) +1;
			#$exon = 'ATCTTGGGTCTTT';
			while($exon =~ /(?=$g)/g){
				push @pos, pos $exon;
				$total_match +=1;
				for my $q(1..length($g)){
					$matches[$pos[$#pos]+$q-1] +=1;
				}
			}
		}
		## deal with the matching positions
		## separate the regions by removing gaps larger than 6bp
		my $match_string = join('',@matches);
		$match_string =~ s/^0+//g;
		$match_string =~ s/0+$//g;
		my @clusters = split(/[0]{$gap,}/, $match_string);
		my $score = 0;
		
		foreach my $c (@clusters){
			my $single_score = 0;
			my $short_int = 0;
			my $long_int  = 0;
			#$c = '111100011100111000000111100000111';
			$short_int = ()= $c =~ /(?=[1-9]+[0]{1,3}[1-9]+)/g;
			$long_int  = ()= $c =~ /(?=[1-9]+[0]{4,6}[1-9]+)/g;
			
			my @set = split(/[0]+/, $c);
			foreach my $s (@set){
				my @bits = split(//, $s);
				my $max = (sort {$b <=> $a} @bits)[0];
				if ($max >= 1){
					$single_score += 2**($max-1) ;
				}
			}
			
			if ($#set > 0){
				#$score += 2**($#set);
				$single_score += (2*$short_int + $long_int ) + $#set+1;
			}else{
				$single_score += $#set +1;
			}

			$single_score /= ($#group +1);
			$score += $single_score;
		}
		if ($#clusters < 0){
			$score = 0;
		}else{
			$score *= length($group[0]);
			$score = log($score)/log(2);
			$score /= length($exon);
			$score *= 100;
			$score = sprintf("%.3f", $score);
		}
		#$$count{$line[0]}{$m}{'AI1'} = $count_u;
		$$count{$line[0]}{$m}{'A'}   = $score;
		#$$count{$line[0]}{$m}{'AI2'} = $count_d;
	}
}

sub hexamer_score
{
	my $ref   = shift;
	my $mut   = shift;
	my $motif = shift;
	my $esr = shift;
	my $gain_loss = shift;
	my $snp = shift;
	my $snp_id = shift;


	$$ref =~ /^([acgtn]+)([ACGTN]+)([acgtn]+)$/;
	my $ref_exon = substr($1,-8).$2.substr($3,0,8);
	$$mut =~ /^([acgtn]+)([ACGTN]+)([acgtn]+)$/;
	my $mut_exon = substr($1,-8).$2.substr($3,0,8);
	
	$ref_exon = uc($ref_exon);
	$mut_exon = uc($mut_exon);
	my %motif_ref_mut = ();
	my %motif_quant = ();
	
	my $mask = $ref_exon ^ $mut_exon;
	$mask =~ /[^\0]/g;
	my $snp_pos = $-[0];

	my @ref_match_loc = ();
	my @mut_match_loc = ();
	foreach my $m (sort keys %$motif){
		my $affected_region_ref = substr($ref_exon, $snp_pos-length($m)+1,2* length($m)-1);
		my $affected_region_mut = substr($mut_exon, $snp_pos-length($m)+1,2* length($m)-1);
		while($affected_region_ref =~ /(?=$m)/g){
			if ($$motif{$m} ne 'chasin-ess' and $$motif{$m} ne 'chasin-ese' and $$motif{$m} ne 'burge-ess' and $$motif{$m} ne 'burge-ese'){
				$motif_ref_mut{"$$motif{$m}:$-[0]"}{'ref'} = 1; ## has more than one binding seq, so use motif name.
				#$motif_ref_mut{"$m:$-[0]"}{'ref'} = 1;
				$motif_quant{$$motif{$m}}{'ref'} +=1;
				#$motif_quant{$m}{'ref'} +=1;
			}else{## chasin/burge ESE/ESS
				$motif_ref_mut{"$m:$-[0]"}{'ref'} += 1;
				$motif_quant{$m}{'ref'} +=1;
			}			
		}
		while($affected_region_mut =~ /(?=$m)/g){
			#push @mut_match_loc, "$m:$-[0]";
			if ($$motif{$m} ne 'chasin-ess' and $$motif{$m} ne 'chasin-ese' and $$motif{$m} ne 'burge-ess' and $$motif{$m} ne 'burge-ese'){
				$motif_ref_mut{"$$motif{$m}:$-[0]"}{'mut'} += 1;
				#$motif_ref_mut{"$m:$-[0]"}{'mut'} += 1;
				$motif_quant{$$motif{$m}}{'mut'} +=1;
				#$motif_quant{$m}{'mut'} +=1;
			}else{
				$motif_ref_mut{"$m:$-[0]"}{'mut'} += 1;
				$motif_quant{$m}{'mut'} +=1;
			}
		}
	}
	foreach my $m (keys %motif_ref_mut){## only used for motif-affected-location. Not for the count change purpose
		if ($motif_ref_mut{$m}{'ref'} eq 1 and $motif_ref_mut{$m}{'mut'} eq ''){
			#Loss
			my @info = split(/:/, $m);
			#my $off = $snp_pos - $info[1]+1;
			my $off = -1;
			if ($info[0] =~ /^[ACGTN]+/){
				$off = length($info[0]) - $info[1];
			}else{
				my @group = split(/\|/, $known_motif{$info[0]});
				$off = length($group[0]) -$info[1];				
			}
			
			$$esr{'loss'}{$info[0]}{$off} +=1;
		}elsif ($motif_ref_mut{$m}{'ref'} eq '' and $motif_ref_mut{$m}{'mut'} eq 1){
			#Gain
			my @info = split(/:/, $m);
			my $off = -1;
			if ($info[0] =~ /^[ACGTN]+/){
				$off = length($info[0]) - $info[1];
			}else{
				my @group = split(/\|/, $known_motif{$info[0]});
				$off = length($group[0]) -$info[1];
			}
			$$esr{'gain'}{$info[0]}{$off} +=1;
		}else{}
	}
	
	foreach my $m (keys %motif_quant){## used for motif-count-change. Record the count change in the SNP affected region. 
		if ($motif_quant{$m}{'ref'} eq ''){
			$motif_quant{$m}{'ref'} =0;
		}
		if ($motif_quant{$m}{'mut'} eq ''){
			$motif_quant{$m}{'mut'} = 0;
		}
	}
	
	$$snp{$snp_id}{'loss'}=0;
	$$snp{$snp_id}{'gain'}=0;
	foreach my $m (keys %motif_quant){
		#my $typ = '';
		#if ($m =~ /^[ACGTN]+/){
		#	if ($$motif{$m} =~ /ess|ese/){
		#		$typ = $$motif{$m};
		#		$typ =~ s/^chasin-|burge-//;
		#	}else{
		#		$typ = $$motif{$m};
		#	}
		#}else{
		#	#my @mo = grep {$$motif{$_} eq $m} keys %$motif;
		#	$typ = $m;
		#}
		if ($motif_quant{$m}{'ref'} > $motif_quant{$m}{'mut'}){
			#$$gain_loss{"loss:$typ"}{$m} +=1;
			$$gain_loss{"loss:$m"} += $motif_quant{$m}{'ref'} - $motif_quant{$m}{'mut'} ;
			$$snp{$snp_id}{'loss'} +=1;
		}elsif($motif_quant{$m}{'ref'} < $motif_quant{$m}{'mut'}){
			#$$gain_loss{"gain:$typ"}{$m} +=1;
			$$gain_loss{"gain:$m"} += $motif_quant{$m}{'mut'} - $motif_quant{$m}{'ref'} ;
			$$snp{$snp_id}{'gain'} +=1;
		}else{}
	}
}
