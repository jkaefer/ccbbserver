#!/usr/bin/perl
use warnings;
use strict;

my ($uploaded_fl, $input_desp, $ts, $flhand ) = @ARGV;
open(RE, ">>$flhand")||die "$!";

#print RE "perl ./source/prepare_features.pl $uploaded_fl $input_desp ./query_jobs/job_$ts";

system("perl ./source/prepare_features.pl $uploaded_fl $input_desp ./query_jobs/job_$ts");
#print RE "   success.1       $?  ";


##----------------------------------------- output prediction result --------------------------------------------
my $var_c_eb = `wc -l < ./query_jobs/job_$ts/$input_desp.weka.data.inner.info`;
my $var_c_ss = `wc -l < ./query_jobs/job_$ts/$input_desp.weka.data.ss.info`;
chomp($var_c_eb); chomp($var_c_ss);
my $eb_download_link = "#";
my $ss_download_link = "#";

#print RE "</br>$var_c_eb\n$var_c_ss\n";


if($var_c_eb >0){
	$eb_download_link = "http://watson.compbio.iupui.edu/regSNP-splicing/query_jobs/job_$ts/$input_desp.exonic_snp_disease_causing_probabilty.exonBody";
}
if($var_c_ss >0){
    $ss_download_link = "http://watson.compbio.iupui.edu/regSNP-splicing/query_jobs/job_$ts/$input_desp.exonic_snp_disease_causing_probabilty.splicingSite";
}

print RE "<div><p><mark><b>Download your prediction result here:</b></mark><p>";
print RE"<p><a href= $eb_download_link>    #SNVs on exon body:      $var_c_eb</a></p>\n";
print RE"<p><a href= $ss_download_link>    #SNVs on splicing sites: $var_c_ss</a></p></div><br>\n";



#$var_c_ss =1;
#$input_desp = 'clinvar.pathogenic';
#$uploaded_fl = './query_jobs/job_2015-4-22-17-0-49/clinvar.pathogenic';
#$ts = '2015-4-22-17-0-49';


## ------------------------------------------------------ output prediction result for SNV on splicing site -------------------------------------
if($var_c_ss >0){
    my $status = system ("perl ./source/code_package/disease_causing_prob_separate.pl $input_desp ss $uploaded_fl ./query_jobs/job_$ts");
#	print RE "perl ./source/code_package/disease_causing_prob_separate.pl $input_desp ss $uploaded_fl ./query_jobs/job_$ts";
#	print RE "<p>$status----------------------------------------------</p>";
#	print RE "<br/><br/><br><hr><br>";
	print RE "<table id=\"splicingSite\"> \n<caption>Predictions for input SNVs on <b><mark>Splicing sites</mark></b></br>Disease causing probability Cut-off:</br><font style=\"background-color:#d0e3f0\">0.44 (Pathogenic Possible: FDR = 0.15, Sensitivity=0.91, Specifity=0.82,MCC=0.73 )</font><br><font style= \"background-color:#ff99b3\"> 0.67 (Pathogenic Probable: FDR = 0.1, Sensitivity=0.82, Specifity=0.89, MCC=0.7)</font></caption>";

	open (OUT, "./query_jobs/job_$ts/$input_desp.exonic_snp_disease_causing_probabilty.splicingSite")||die "$!";
	my $h1 = <OUT>;
	$h1 =~ s/\s+$//;
	my @head = split(/\t|\s+/, $h1);
	print RE "<tr>";
	foreach my $h (0..$#head){
	    print RE "<th>$head[$h]</th>";
	}
	print RE "</tr>";
	while(<OUT>){
	        s/\s+$//;
	        my @line = split(/\t|\s+/, $_);
	        print RE "<tr>";
	        foreach my $e(0..$#line){
	                if($line[4] >= 0.67 ){ # probable, FDR
	                    print RE "<td id =\"probable\">$line[$e]</td>";
	                }elsif($line[4] >= 0.44 ){ ## possible, FDR < 
			    print RE "<td id = \"possible\">$line[$e]</td>";
			}else{
	                    print RE "<td>$line[$e]</td>";
	                }
	        }
	        print RE "</tr>";
	}
	
	print RE "</table>";
	close OUT;
}
##----------------------------------------- output prediction result for SNV on exonBody--------------------------------------

if($var_c_eb >0){
	my $status = system("perl ./source/code_package/disease_causing_prob_separate.pl $input_desp inner $uploaded_fl ./query_jobs/job_$ts");
#	print RE "<p>perl ./source/code_package/disease_causing_prob_separate.pl $input_desp inner $uploaded_fl ./query_jobs/job_$ts 2>&1</p>";
#	print RE "<p>$status</p>";
	print RE "<br><hr><br>";
	print RE "<table id=\"exonBody\"> \n<caption>Predictions for input SNVs on <mark><b>exon body</b></mark></br>Disease causing probability Cut-off:</br><font style=\"background-color:#d0e3f0\">0.48 (Patheogenic Possible: FDR=0.25, Sensitivity = 0.76, specifity = 0.76,MCC=0.52)</font><br><font style= \"background-color:#ff99b3\"> 0.6 (Pathogenic Probable: FDR = 0.1, Sensitivity = 0.53, Specifity = 0.88, MCC=0.43)</font></caption>";

	open (OUT, "./query_jobs/job_$ts/$input_desp.exonic_snp_disease_causing_probabilty.exonBody")||die "$!";
	my $h1 = <OUT>;
	$h1 =~ s/\s+$//;
	my @head = split(/\t|\s+/, $h1);
	print RE "<tr>";
	foreach my $h (0..$#head){
	        print RE "<th>$head[$h]</th>";
	}
	print RE "</tr>";
	while(<OUT>){
	        s/\s+$//;
	        my @line = split(/\t|\s+/, $_);
	        print RE "<tr>";
	        foreach my $e(0..$#line){	
	                if($line[4] >= 0.6 ){ # probable, FDR
	                    print RE "<td id =\"probable\">$line[$e]</td>";
	                }elsif($line[4] >= 0.48){ ## possible, FDR < 
			    print RE "<td id = \"possible\">$line[$e]</td>";
			}else{
	                    print RE "<td>$line[$e]</td>";
	                }
	        }
	        print RE "</tr>";
	}	
	print RE "</table>";
	close OUT;
}


close OUT;


print RE "<br>";
#print RE "<p>***</p>";
print RE "<br />";

print RE "</body>";
close RE;

rename("./query_jobs/job_$ts/result.html", "./query_jobs/job_$ts/result.tmp")||die "$!";
rename("./query_jobs/job_$ts/result", "./query_jobs/job_$ts/result.html")||die "$!";
