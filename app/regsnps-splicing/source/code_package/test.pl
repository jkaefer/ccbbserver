#!/usr/bin/perl
use warnings;
use strict;


my @abc = (5,7,6,2,7,1);
my @sorted_list = sort {$abc[$a] <=> $abc[$b]} 0..$#abc;
my @idx = @abc[@sorted_list];

my $job = '';
my $desp = '';
my $type = '';
my %snp = ();
my @available_snp = ();
my $snpid = '';

#print O1 "3 done\n";
my @ss_stat = ();
my @eb_stat = ();
open(SS, "./regsnps-splicing/model/stat.ss")||die "$!";
while(<SS>){
        s/\s+$//;
        my @line = split(/,/, $_);
        push @ss_stat, [@line];
}
close SS;

open(EB, "./regsnps-splicing/model/stat.exonbody")||die "$!";
while(<EB>){
        s/\s+$//;
        my @line = split(/,/, $_);
        push @eb_stat, [@line];
}
close EB;

$type = 'exonBody';
#open(PRB, "$job/$desp.$type.prediction")||die "$!";
#while(<PRB>){
#       s/\s+$//;
#        s/^\s+//;
#        if(/^\d/){
#                my @line = ();
#                if(/\+/){
#                        s/\+//;
#                        @line = split(/\s+/, $_);
#                }else{
#                        @line = split(/\s+/, $_);
#                        $line[$#line] = 1- $line[$#line];
#                }
#                $snp{$available_snp[$snpid]}{'prob'} = $line[$#line];
                #$snpid +=1;
                my $prob = 0;
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
                        my $diff = abs($prob - $cutoff[$#cutoff]);
                        if ($diff <= $mindiff){
                                $mindiff = $diff;
                                $mindi_idx = $i
                        }
                }
                @cutoff = @{$stat[$mindi_idx]};
                my $fdr = $cutoff[2]/($cutoff[0]+$cutoff[2]);$fdr = sprintf("%.3f",$fdr);
                #$snp{$available_snp[$snpid]}{'FDR'} = $fdr;
                #$snpid +=1;
#       }
#}
#close PRB;
		print "$fdr\t$mindi_idx\t$mindiff\n";
