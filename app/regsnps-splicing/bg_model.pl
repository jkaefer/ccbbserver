#!/usr/bin/perl
use warnings;
use strict;

print "PROBLEM bg \n";

my ($uploaded_fl, $input_desp, $ts, $flhand ) = @ARGV;
open(RE, ">>$flhand")||die "$!";

#print RE "perl ./regsnps-splicing/source/prepare_features.pl $uploaded_fl $input_desp ./regsnps-splicing/query_jobs/job_$ts";
my $output=qx(pwd);
print $output;
system("perl ./regsnps-splicing/source/prepare_features.pl $uploaded_fl $input_desp ./allJobs/Splicejob_$ts");
#print RE "   success.1       $?  ";


##----------------------------------------- output prediction result --------------------------------------------
my $var_c_eb = `wc -l < ./allJobs/Splicejob_$ts/$input_desp.weka.data.inner.info`;
my $var_c_ss = `wc -l < ./allJobs/Splicejob_$ts/$input_desp.weka.data.ss.info`;
chomp($var_c_eb); chomp($var_c_ss);
my $eb_download_link = "#";
my $ss_download_link = "#";

#print RE "</br>$var_c_eb\n$var_c_ss\n";


if($var_c_eb >0){
	$eb_download_link = "./allJobs/Splicejob_$ts/$input_desp.exonic_snp_disease_causing_probabilty.exonBody";
}
if($var_c_ss >0){
    $ss_download_link = "./allJobs/Splicejob_$ts/$input_desp.exonic_snp_disease_causing_probabilty.splicingSite";
}

my $all_download_link = "./allJobs/Splicejob_$ts/$input_desp.exonic_snp_disease_causing_probabilty.combine";
my $invalid_var_download_link = "./allJobs/Splicejob_$ts/$input_desp.invalid_input";
print RE "<div><p><mark><b>Download your prediction result & analysis log here:</b></mark><p>";
#print RE"<p><a href= $eb_download_link>    #SNVs on exon body:      $var_c_eb</a></p>\n";
#print RE"<p><a href= $ss_download_link>    #SNVs on splicing sites: $var_c_ss</a></p></div><br>\n";
print RE"<p><a href= $all_download_link>Download prediction result</a>,  <a href=$invalid_var_download_link > View invalid input variant</a> </p></div>\n";
#print RE"<p><a href= $invalid_var_download_link>View invalid input variant</a></p></div><br>\n";
print RE "<p>We suggest using FDR as cut-off: FDR<0.05 as deleterious, 0.05<=FDR<0.5 as potential deleterious, and 0.5<=FDR<=1 as benign.</p>\n";



## ------------------------------------------------------ output prediction result for SNV on splicing site -------------------------------------
my %combined_pred = ();
my @snp_order = ();
open(IN,"$uploaded_fl" )||die "$!"; ## 'query' or uploaded file
while(<IN>){
	s/\s+$//;
	my @line = split(/\s+/, $_);
	#$snp{$line[0]}{'type'} = $line[3];
	if ($#line == 2){
		push @snp_order, "$line[0]:$line[1]-$line[2]";
	}elsif($#line > 2 and $line[3] =~ /[ACGT]+/ and $line[4] =~ /[ACGT]+/){
		my @alta = split(/,/, $line[4]);
		foreach my $al (@alta){
			push @snp_order, "chr$line[0]:$line[1]:$line[3]-$al"; 
		}
	}
}
close IN;

my $fl1 = "./allJobs/Splicejob_$ts/$input_desp.exonic_snp_disease_causing_probabilty.splicingSite";
my $fl2 = "./allJobs/Splicejob_$ts/$input_desp.exonic_snp_disease_causing_probabilty.exonBody";
my $fl3 = "./allJobs/Splicejob_$ts/$input_desp.exonic_snp_disease_causing_probabilty.combine";
my $h = '';
#print RE "ok 1\n";
if($var_c_ss >0){
    print "PROBLEM sep \n";
    my $status = system ("perl ./regsnps-splicing/source/code_package/disease_causing_prob_separate.pl $input_desp ss $uploaded_fl ./allJobs/Splicejob_$ts");
	#print RE "perl ./regsnps-splicing/source/code_package/disease_causing_prob_separate.pl $input_desp ss $uploaded_fl ./regsnps-splicing/query_jobs/job_$ts\n";
	#print RE "$status\n";
	
	open(PD, $fl1)||die "$!";
	#print RE "open ok\n";
	$h = <PD>;
	while(<PD>){
		s/\s+$//;
		#my @line = split(/\s+/, $_);
		my @line = split(',', $_);
		$combined_pred{"$line[0]:$line[1]:$line[2]-$line[3]"} = $_;
	}
	close PD;

}
#print RE "ok 2\n";
if($var_c_eb >0){
	print "PROBLEM combine \n";
	my $status = system("perl ./regsnps-splicing/source/code_package/disease_causing_prob_separate.pl $input_desp inner $uploaded_fl ./allJobs/Splicejob_$ts");
	
	#print RE "$status\n";
	open(PD, $fl2)||die "$!";
	$h = <PD>;
	while(<PD>){
		s/\s+$//;
		#my @line = split(/\s+/, $_);
		my @line = split(',', $_);
		if(!/NA/){
		$combined_pred{"$line[0]:$line[1]:$line[2]-$line[3]"} = $_;
		}
	}
	close PD;
}

my @adj_fdr  =();
my @disp = ();
my %f_idx = ();
open(AJ, ">$fl3.adj")|| die "$!";
open(OT, ">$fl3")||die "$!";
print OT "$h";
print AJ "$h";
foreach my $so (@snp_order){
	if($combined_pred{$so} ne ''){
		print OT  "$combined_pred{$so}\n";
                my @data = split(/,/,$combined_pred{$so});
                push @adj_fdr, $data[5];
                push @disp, $data[4];

	}
}

close OT;
my @f_idx_map = ();
@adj_fdr = sort {$a <=> $b} @adj_fdr;
@f_idx_map = sort {$disp[$b] <=> $disp[$a]} 0..$#disp;
#@f_idx = reverse @f_idx;
foreach my $f (0..$#f_idx_map){
        $f_idx{$f_idx_map[$f]} = $f;
}

open(OT, $fl3)||die "$!";
my $tmpd = <OT>;
my $lc = 0;
#print AJ " @disp\t@f_idx\t@adj_fdr\t@adj_fdr[@f_idx]\n";
while(<OT>){
	s/\s+$//;
        my @data = split(/,/, $_);
        print AJ "$data[0],$data[1],$data[2],$data[3],$data[4],";
        #print AJ "$adj_fdr[($#f_idx) - $f_idx[$lc]]";
	print AJ "$adj_fdr[ $f_idx{$lc}]";
        for my $i (6..$#data){
                print AJ ",$data[$i]";
        }
	print AJ "\n";
	$lc +=1;
}
close OT;
close AJ;

#convert to json for javascript tables
system("python ./regsnps-splicing/conv.py $fl3 ./allJobs/Splicejob_$ts/modded.$input_desp.exonic_snp_disease_causing_probabilty.combine");


	#print RE "<table id=\"splicingSite\"> \n<caption>Predictions for input SNVs on <b><mark>Splicing sites</mark></b></br>Disease causing probability Cut-off:</br><font style=\"background-color:#d0e3f0\">0.44 (Pathogenic Possible: FDR = 0.15, Sensitivity=0.91, Specifity=0.82,MCC=0.73 )</font><br><font style= \"background-color:#ff99b3\"> 0.67 (Pathogenic Probable: FDR = 0.1, Sensitivity=0.82, Specifity=0.89, MCC=0.7)</font></caption>";







