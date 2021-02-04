#!/usr/bin/perl
use warnings;

print "PROBLEM disease \n";
  
my $desp  = $ARGV[0];
my $type  = $ARGV[1];
my $qfile = $ARGV[2];
my $job   = $ARGV[3];

my $infile = "$job/$desp.snp.ref.mut.seq.proximity";
my %snp = ();

#open(O1, ">s1")||die "$!";
#print O1 "0 done\n";
my @snp_order = ();
open(IN,"$qfile" )||die "$!"; ## 'query' or uploaded file
while(<IN>){
	s/\s+$//;
	s/^chr//;
	my @line = split(/\s+/, $_);
	#$snp{$line[0]}{'type'} = $line[3];
	if ($#line == 2){
		$snp{"$line[0]:$line[1]-$line[2]"}{'prob'} = -1; ## all initialized to be -1. Some will not have prediction
		$snp{"$line[0]:$line[1]-$line[2]"}{'ss'}   = -1; ## exonbody or splicing site
		push @snp_order, "$line[0]:$line[1]-$line[2]";
	}elsif($#line > 2 and $line[3] =~ /[ACGT]+/ and $line[4] =~ /[ACGT]+/){
		my @alta = split(/,/, $line[4]);
		foreach my $al (@alta){
		    $snp{"chr$line[0]:$line[1]:$line[3]-$al"}{'prob'} = -1;
		    $snp{"chr$line[0]:$line[1]:$line[3]-$al"}{'ss'}   = -1;
                    push @snp_order, "chr$line[0]:$line[1]:$line[3]-$al"; 
		}
	}
}
close IN;


#print O1 "1 done\n";
my %flist_ss = ();
my %flist_eb = ();
open(ON, "./regsnps-splicing/source/code_package/model/weka.head.ss")||die "$!";
my $lc = 0;
while(<ON>){
    s/\s+$//;
    if(/real$/){
	my @line = split(/\s+/, $_);
	$flist_ss{$lc} = $line[1];
	$lc +=1;
    }
}
close ON;

$lc = 0;
open(ON, "./regsnps-splicing/source/code_package/model/weka.head.inner")||die "$!";
while(<ON>){
    s/\s+$//;
    if(/real$/){
	my @line = split(/\s+/, $_);
	$flist_eb{$lc} = $line[1];
	$lc +=1;
    }
}
close ON;


open(IN, $infile) ||die "$!";
while(<IN>){
	s/\s+$//;
	 my @line =  split(/\t/, $_);
	 $snp{$line[2]}{'transcript'} = $line[1];
	#$snp{$line[2]}{'maf'} = $line[$#line];
}
close IN;

my $offset = 0;
if($type eq 'ss'){
	$offset = 14;
}else{
	$offset = 15;
}
my @available_snp = ();

#print O1 "2 done\n";

## +++++++++++++++++++++ check which SNPs have been predicted  ++++++++++++++++++++++++
open(IN, "$job/$desp.weka.data.$type.info")||die "$!";
while(<IN>){
	s/\s+$//;
	my @line = split(/,/, $_);
	push @available_snp, $line[$offset];
	$snp{$line[$offset]}{'ss'} = 1;
	if($type eq 'ss'){
	    for my $i (0..$#line){
		if($i < $offset){
                    $snp{$line[$offset]}{$flist_ss{$i}} = $line[$i];
		}elsif($i > ($offset +3 ) ){
		    $snp{$line[$offset]}{$flist_ss{$i-4}} = sprintf("%.2f", $line[$i]);
		}else{}
	    }
	    $snp{$line[$offset]}{'proximity_3'} = $line[15];
	    $snp{$line[$offset]}{'proximity_5'} = $line[16];
	    $snp{$line[$offset]}{'loc'} = "ON";
	}else{ ## exon body		`
	    for my $i (0..$#line){
		if($i < 5){
		    $snp{$line[$offset]}{$flist_eb{$i}} = $line[$i];
		}elsif($i >= 7 and $i < 15){
		    $snp{$line[$offset]}{$flist_eb{$i-2}} = $line[$i];
		}elsif($i == 16 or $i ==17){
		    $snp{$line[$offset]}{$flist_eb{$i-3}} = $line[$i];
		}elsif($i >= 19){
		    $snp{$line[$offset]}{$flist_eb{$i-4}} = sprintf("%.2f",$line[$i]);
		}else{}
	    }
	    $snp{$line[$offset]}{'junction_score_change'} = 0;
	    $snp{$line[$offset]}{'loc'} = "OFF";
	}
}
close IN;

my $snpid = 0;
if($type eq 'inner'){
	$type = 'exonBody';
}
if ($type eq 'ss'){
	$type = "splicingSite";
}

#print O1 "3 done\n";
my @ss_stat = ();
my @eb_stat = ();
open(SS, "./regsnps-splicing/source/code_package/model/stat.4.ss")||die "$!";
while(<SS>){
	s/\s+$//;
	my @line = split(/,/, $_);
	push @ss_stat, [@line];
}
close SS;

open(EB, "./regsnps-splicing/source/code_package/model/stat.5.exonbody")||die "$!";
while(<EB>){
	s/\s+$//;
	my @line = split(/,/, $_);
	push @eb_stat, [@line];
}
close EB;

my @adjust_fdr = ();
open(PRB, "$job/$desp.$type.prediction")||die "$!";
while(<PRB>){
	s/\s+$//;
	s/^\s+//;
	if(/^\d/){
		my @line = ();
		if(/\+/){
			s/\+//;
			@line = split(/\s+/, $_);
		}else{
			@line = split(/\s+/, $_);
			$line[$#line] = 1- $line[$#line];
		}
		$snp{$available_snp[$snpid]}{'prob'} = $line[$#line];
		#$snpid +=1;
		
		my @stat = ();
		if($type eq 'exonBody'){
			@stat = @eb_stat;
		}else{
			@stat = @ss_stat;
		}

		my @fdr_diff = ();
		my $mindiff = 100;
		my $mindi_idx = -1;
		my @cutoff = ();
		for my $i (0..$#stat){
			@cutoff = @{$stat[$i]};
			my $diff = abs($snp{$available_snp[$snpid]}{'prob'} - $cutoff[$#cutoff]);
			if ($diff <= $mindiff){
				$mindiff = $diff;
				$mindi_idx = $i
			}
		}
		@cutoff = @{$stat[$mindi_idx]};
		#my $fdr = $cutoff[2]/($cutoff[0]+$cutoff[2]);$fdr = sprintf("%.3f",$fdr);
		my $fdr = ((1+$cutoff[0])/(2+$cutoff[0]+$cutoff[1]))/((1+$cutoff[2])/(2+$cutoff[2]+$cutoff[3]));$fdr = sprintf("%.3f",1/$fdr);
		$snp{$available_snp[$snpid]}{'FDR'} = $fdr;
		push @adjust_fdr, $fdr;
		$snpid +=1;
	}
}
close PRB;


#print O1 "4 done\n";

my @user_att = ();
my @user_common_name = ();
open(CAT, "$job/checked_attributes")||die "$!";
while(<CAT>){
	s/\s+$//;
if(/exonIntron/){
	push @user_att, "exon_intron_len_1","exon_intron_len_2","exon_intron_len_3";
	push @user_common_name, "upstream intron length","exon length","downstream intron length";
}elsif(/proximity/){
	push  @user_att, "proximity_3", "proximity_5";
	push @user_common_name, "proximity acceptor site","proximity donor site";
}elsif(/evolution/){
	push  @user_att, "evolution";
	push @user_common_name, "conservation";
}elsif(/junction/){
	push  @user_att, "junction_score_1", "junction_score_2";
	push @user_common_name, "acceptor_site strength","donor_site strength";
}elsif(/junc_change/){
	push  @user_att, "junction_score_change";
	push @user_common_name, "junction score change";
}elsif(/sfrs/){
	push  @user_att, "sfrs1_3", "sfrs1_5", "sfrs2_3", "sfrs2_5", "sfrs5_3", "sfrs5_5", "sfrs6_3", "sfrs6_5";
	push @user_common_name, "SFRS1(ref)","SFRS1(mut)","SFRS2(ref)","SFRS2(mut)","SFRS5(ref)","SFRS5(mut)","SFRS6(ref)", "SFRS6(mut)";
}elsif(/cluster_score/){
	push  @user_att, "cluster_score_change", "cluster_score_ref"; 
	push @user_common_name, "cluster score(change)","cluster score(original)";
}elsif(/ef/){
	push  @user_att, "EF_score";
	push @user_common_name, "EF score"
}elsif(/rna2nd/){
	push  @user_att, "RNA_2nd_struct";
	push @user_common_name, "RNA_2nd structure (change)";
}elsif(/PTM/){
	push  @user_att, "ptm_1";
	push @user_common_name, "PTM";
}elsif(/Pfam/){
	push  @user_att, "pfam_coverage";
	push @user_common_name, "Pfam";
}elsif(/disorder/){
	push  @user_att, "disorder_score_12";
	push @user_common_name, "ave. disorder score";
}elsif(/ss/){
	push  @user_att, "ss_3";
	push @user_common_name, "ave. structure probability";
}elsif(/asa/){
	push  @user_att, "asa_1";
	push @user_common_name, "ave. ASA";
}else{}
}
close CAT;


open(OUT, ">$job/$desp.exonic_snp_disease_causing_probabilty.$type")||die "$!";
print OUT "chr,pos,ref,alt,disease prob,FDR,transcript,ON/OFF splicing site";
# exon_length proximity_2_acceptor proximity_2_donor junction_score_acceptor junction_score_donor junction_score_change evolution ave_disorder_score ave_ASA ave_SS_coil\n";
#foreach my $att (@user_att){
foreach my $cm(@user_common_name){
	print OUT ",$cm";
}
print OUT "\n";

my @disp = ();
foreach my $ss (keys %snp){
	push @disp, $snp{$ss}{'prob'};
}

my @fdr_rank = sort {$disp[$b] <=> $disp[$a]} 0..$#disp;

my $f_idx = 0;
foreach my $ss (keys %snp){
        my @loc = split(/:|-/, $ss);
        print "$snp{$ss}{'prob'}";
	if($snp{$ss}{'prob'} != -1){
		my $a = $loc[1]-1;
		#print OUT "$loc[0] $loc[1] $loc[2] $loc[3] $snp{$ss}{'prob'} $snp{$ss}{'FDR'} $snp{$ss}{'transcript'} $snp{$ss}{'exon_length'} $snp{$ss}{'proximity_2_acceptor'} $snp{$ss}{proximity_2_donor} $snp{$ss}{'junction_score_acceptor'} $snp{$ss}{'junction_score_donor'} $snp{$ss}{'junction_score_change'} $snp{$ss}{'evolution'} $snp{$ss}{'ave_disorder_score'} $snp{$ss}{'ave_ASA'} $snp{$ss}{'ave_SS_coil'}\n";
		print OUT "$loc[0],$loc[1],$loc[2],$loc[3],$snp{$ss}{'prob'},$snp{$ss}{'FDR'},$snp{$ss}{'transcript'},$snp{$ss}{'loc'}";
		foreach my $att (@user_att){
			print OUT ",$snp{$ss}{$att}";
		}
		print OUT "\n";
	}else{
		print OUT "$loc[0],$loc[1],$loc[2],$loc[3],NA,NA,$snp{$ss}{'transcript'}";
		foreach my $att (@user_att){
		    print OUT ",NA";
		}
		print OUT "\n";
	}
	$f_idx +=1;
}

close OUT;
#close MAF;
#print O1 "5 done\n";

#close O1;