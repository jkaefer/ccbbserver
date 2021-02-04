#!usr/bin/perl -w

# calculate logor and maxOR for each SNP

my $job = $ARGV[1];
my $motiflist = "./source/code_package/db/pwm.all.list";
my @pwmrows = ("A","C","G","T");
my %motifs = (); # for all pwm

#my $infile = "neutral.snp.ref.mut.seq";
my $infile = $ARGV[0];
open (OUT ,">$job/$infile.motif_change.out")||die "$!";
print OUT "PWM\tid\tchr\tstrand\tindel_left\tindel_right\tgene\tref\talt\tindel\tMut.score\tMut.pos\tRef.score\tRef.pos\n";

#read in all pssm
open (IN,"$motiflist")||die "$!";
while(<IN>)
{
   s/\s+$//;
    open (PW,"$_")||die "$!";
    my $flname = $_;
    
    my $cc  = 0;
    my %pwm = (); # single pssm
    my $mo_len = 0;
    while(<PW>)
    {
		s/\s+$//;
		s/^\s+//;
		my @line = split(/\s+/,$_);
		$mo_len = $#line +1;
		$pwm{$pwmrows[$cc]} = [@line];
		$cc ++;
   }
   $pwm{'len'} = $mo_len;
   close PW;
   $flname =~ s/source\/code_package\/db\/motif\/PWM\///;
   $motifs{$flname} = {%pwm};
   %pwm = ();
    
}
close IN;

my @pwm_out_ord = ();
open(IN, "./source/code_package/db/pwm_output_order")||die "$!";
while(<IN>){
s/\s+$//;
push @pwm_out_ord, $_;
}

close IN;


open(IN, "$job/$infile")||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	#$line[3] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	#my $ref_exon = $2;
	#$line[4] =~ /([acgtn]+)([ACGTN]+)([acgtn]+)/;
	#my $mut_exon = $2;
	
	$ref_exon = $line[3];
	$mut_exon = $line[4];
	$ref_exon = uc($ref_exon);
	$mut_exon = uc($mut_exon);
	
	my $mask = $ref_exon ^ $mut_exon;
	while($mask =~ /[^\x00]+/g)
	{
		my @fq = ();
		$fq[0] = "$line[0]\t$line[1]\t$line[2]\t".substr($ref_exon, $-[0],1)."\t".substr($mut_exon, $-[0],1)."\tINDEL";
		$fq[1] =  substr($ref_exon, 0, $-[0])."-".substr($mut_exon, $-[0],1)."-".substr($ref_exon, $-[0]+1);
		searchReg(\@fq,\%motifs,\%motif_change);
		
		#foreach my $pp (keys %motif_change)
		foreach  my $pp (@pwm_out_ord)
		{
		    my $m_pp = $pp;
		    #$m_pp =~ s/\/data\/pro\/as\/data\/PWMDir\///;
		    #$m_pp =~ s/\/data\/pro\/as\/data\/pwms_all_motifs\///;
		    $m_pp =~ s/source\/code_package\/db\/motif\/PWM\///;
		    
		    foreach my $bp (keys %{$motif_change{$pp}})
		    {
				#print OUT "$pp\t$bp\t$data_seq{$bp}";
				print OUT "$m_pp\t$bp";
				for my $con ( @{$motif_change{$pp}{$bp}{'mut'}} )
				{
					print OUT "\t$con";
				}
				for my $con ( @{$motif_change{$pp}{$bp}{'ref'}} )
				{
					print OUT "\t$con";
				}	
				print OUT "\n";
		    }
		    
		}
		
		%motif_change = ();
		
	}
}

close IN;

#---------------------------------------------------------------------------#
sub searchReg
{
    my $seq = shift;
    my $motif = shift;
    my $chan = shift;
    
    my @id = split(/\t/, $$seq[0]);
    my $del_seq = "";
    my $indel_del = "";
    my $indel_ins = "";
    my $mut_pat = "";
    
#    if ($id[$#id] =~ /del/ and ! ($id[8] =~ /ins/) ) # del
#    {
#	my @tmp = split(/del/, $id[8]);
#	$del_seq = $tmp[$#tmp];
#	$mut_pat = 'del';
#   }
#    elsif ($id[8] =~ /INDEL/) # indel  ## MODIFICATION for 1000 genome data
#    {
	#$id[4] =~ /del([ATGC]+)ins([ATGC]+)/;
	$indel_del = $id[3];# ref allele
	$indel_ins = $id[4];# mut allele
#	if($id[2] eq "-")
#	{
#	    $indel_del = revcomp(\$indel_del);
#	    $indel_ins = revcomp(\$indel_ins);
#	}
#	$mut_pat = "indel";
#   }
#    elsif ($id[8] =~ /dup/)						#dup
#    {
#	$mut_pat = 'dup';
#    }
#    elsif ($id[8] =~ /ins[AGTC]+/ and !($id[8] =~ /del/))			#ins
#    {
#	$mut_pat = "ins";
#    }
#    else{}
    
    my @seq_bkp = split(/-/,$$seq[1]); # sequence around breakpoint

    foreach my $mo (keys %$motif)
    {
		my $motif_len = $$motif{$mo}{'len'};
		my @up_seq = split(//,$seq_bkp[0]);
		my @dn_seq = split(//,$seq_bkp[2]);
		
		my $up_sch_seq 	= ""; 	# sequence pulled out for searching motif
		my $dn_sch_seq 	= "";
		my $mut_sch_seq	= ""; 	# whole searching sequence, mut
		my $ref_sch_seq	= "";	# whole searching sequence, ref
	
		if ($#up_seq >= $motif_len and $#dn_seq >= $motif_len)
		{
			$up_sch_seq = join('',@up_seq[$#up_seq-($motif_len-2)..$#up_seq]);
			$dn_sch_seq = join('',@dn_seq[0..$motif_len-2]);
			
			#if ($mut_pat eq 'del')
			#{
			#$mut_sch_seq = $up_sch_seq.$dn_sch_seq;
			#$ref_sch_seq = $up_sch_seq.$del_seq.$dn_sch_seq; # ref seq: plugin deleted seq
			#}
			#elsif ($mut_pat eq 'ins' or $mut_pat eq 'dup')
			#{
			#$mut_sch_seq = $up_sch_seq.$seq_bkp[1].$dn_sch_seq; # mut seq: plugin inserted/duplicated seq
			#$ref_sch_seq = $up_sch_seq.$dn_sch_seq;
			#}
			#elsif ($mut_pat eq 'indel')
			#{
			$mut_sch_seq = $up_sch_seq.$seq_bkp[1].$dn_sch_seq; # mut seq: plugin inserted seq
			$ref_sch_seq = $up_sch_seq.$indel_del.$dn_sch_seq;	#ref seq: plugin deleted seq
			#}
			#else{}
		}
		else
		{
			print "err: flanking region too short!\n$$seq[0]\n\n";
		}
		
		my %search = ();
		$search{'mut'} = $mut_sch_seq;
		$search{'ref'} = $ref_sch_seq;
		motifMatch($$seq[0], $motif, $mo, \%search, $chan);
    }
}
    
sub motifMatch
{
	my $id = shift;
	my $motif = shift;
	my $mo = shift;
	my $s_seq = shift; # searching sequences
	my $mo_change = shift; # the motif changes caused by SV
	
	my $mut_sch_seq = "";
	my $ref_sch_seq = "";
	
	my @order = ('mut','ref');
	
	#search for motif
	my %motif_socre = ();
	#foreach my $mo (keys %$motif)
	#{
	    
	    foreach my $con (@order)
	    {
			my $motif_len =  $$motif{$mo}{'len'};
			my $socre = -65536;
			#my $ref_score = -65536;
			my $tmp_score = 0.0;
			my $pos = 0;
		
			my @sch_seq = split(//,$$s_seq{$con});
			for my $i(0..($#sch_seq-$motif_len+1) )# sliding window
			{
		    	for my $j (0..$motif_len-1 )# for each position
			    {
					$tmp_score += $$motif{$mo}{$sch_seq[$i+$j]}[$j];
			    }
		    	if ($tmp_score > $socre )
			    {
					$socre = $tmp_score;
					$pos = $i-$motif_len+1; # distance away from breakpoint (left most breakpoint)
		    	}
			    $tmp_score = 0.0;
			}
		
			$$mo_change{$mo}{$id}{$con} = [$socre,$pos];
		}
}
	    
	    
	    
sub revcomp
{
	my $s = shift;
	my %rev = ();
	$rev{'A'} = 'T';
    $rev{'T'} = 'A';
    $rev{'G'} = 'C';
    $rev{'C'} = 'G';
	    
    my @seq = split('', $$s);
    my $revs = '';
    foreach my $e (@seq)
    {
		$revs .= $rev{$e};
    }
    return (reverse scalar $revs);
}


