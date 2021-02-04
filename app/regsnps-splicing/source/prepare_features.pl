#!usr/bin/perl
use warnings;
use strict;

print "PROBLEM features \n";

my $infile = $ARGV[0];# input SNP file
my $desp   = $ARGV[1]; # description for input file
my $archive = $ARGV[2];

open(LOG, ">$archive/$desp.log")||die "$!";
print LOG "infile:$infile \ndesption:$desp \narchive location:$archive\n";
#print STDOUT "Hello\n";

### get reference seq for each exon with SNPs
print LOG "......(1/9) extracting exon and neighboring intron sequences......\n";
print LOG "perl ./regsnps-splicing/source/code_package/PrepareSeq/extract_target_sequence.pl $infile $desp $archive";
my $status = system("perl ./regsnps-splicing/source/code_package/PrepareSeq/extract_target_sequence.pl $infile $desp $archive");
print LOG "\n$status\n";

# exon/intron length
print LOG "......(2/9) calculating length of exon and neighboring introns......\n";
$status = system("perl ./regsnps-splicing/source/code_package/ExonIntronLen/exon_intron_length.pl $desp $archive");
print  LOG "<p>$status</p>\n";

# SR_protein matching score
print LOG "......(3/9) calculating Splicing factor binding affinity score......\n";
$status = system("perl ./regsnps-splicing/source/code_package/SFProtein/SR_protein.score.pl /N/slate/jkaefer/splicingdb/motif/SRSF1.pssm $desp $archive");
print  LOG "<p>$status</p>\n";
$status = system("perl ./regsnps-splicing/source/code_package/SFProtein/SR_protein.score.pl /N/slate/jkaefer/splicingdb/motif/SRSF2.pssm $desp $archive");
print  LOG "<p>$status</p>\n";
$status = system("perl ./regsnps-splicing/source/code_package/SFProtein/SR_protein.score.pl /N/slate/jkaefer/splicingdb/motif/SRSF5.pssm $desp $archive");
print  LOG "<p>$status</p>\n";
$status = system("perl ./regsnps-splicing/source/code_package/SFProtein/SR_protein.score.pl /N/slate/jkaefer/splicingdb/motif/SRSF6.pssm $desp $archive");
print  LOG "<p>$status</p>\n";

print "MADE IT HERE feats \n";

# junction score
print LOG "......(4/9) calculating junction score of SNPs affected exons......\n";
$status =  system("perl ./regsnps-splicing/source/code_package/JunctionScore/junction_score.pl $desp $archive");
print  LOG "<p>$status</p>\n";

#evolution
print LOG "......(5/9) calculating the average conservation score(Phylop) for +/- 5bp around SNP......\n";
my $output=qx(pwd);
print $output;
$status =  system("perl /N/slate/jkaefer/splicingdb/phylop_score/acess.phylop.pl $desp.snp.ref.mut.seq.mutBed $archive");
print  LOG "<p>$status</p>\n";

#all RBP
print LOG "......(6/9) calculating the effect of distuptions in RNA binding proteins(RBPs)......\n";
$status =  system("perl ./regsnps-splicing/source/code_package/RBP/cal_RBP_matching_score.pl $desp.snp.ref.mut.seq $archive");
print  LOG "<p>$status</p>\n";

print "Pre R 1 \n";

$status =  system("Rscript ./regsnps-splicing/source/code_package/RBP/magnitude_postProb.r $desp.snp.ref.mut.seq $archive");
print  LOG "<p>$status</p>\n";

print "Post R 1 \n";

##cluster score & ESR-HS score
print LOG "......(7/9) calculating the disruption of ESE/ESS in exons......\n";
$status =  system("perl ./regsnps-splicing/source/code_package/ClusterScore/cluster_score.pl $desp.snp.ref.mut.seq $archive");# generate cluster score
print  LOG "<p>$status</p>\n";

print LOG "......(8/10) RNA 2nd structure and structure distance......\n";
$status =  system("perl ./regsnps-splicing/source/code_package/EFScore/compute_ef_score_around_snp.pl $desp.snp.ref.mut.seq $archive");
print  LOG "<p>$status</p>\n";
$status =  system("perl ./regsnps-splicing/source/code_package/RNAstructure/structure_distance.pl $desp $archive");
print  LOG "<p>$status</p>\n";

##PTM, Pfam, disorder score
print LOG "......(9/10) calculating the protein function level information:Pfam, PTM and protein disorder score......\n";
$status =  system("perl ./regsnps-splicing/source/code_package/PfamPTMDisorder/pfam_ptm_disorderScore.pl $desp.snp.ref.mut.seq $archive");
print  LOG "<p>$status</p>\n";

#ASA,SS feature
print LOG "......(10/10) calculating the protein physical characteristics:Secondary structure(alpha-helix, beta-sheet and coil) and Accessible surface area(ASA)......\n";
$status =  system("Rscript ./regsnps-splicing/source/code_package/ASASS/norm_asa_ss.r $desp $archive");
print "Post R 2 \n";
print LOG "<p>$status</p>\n";

##combine data together
print LOG "......combine features together......\n";
#$status =  `Rscript ./regsnps-splicing/source/code_package/pipeline.selectedFeature.r $desp $archive`;
$status =  system("Rscript ./regsnps-splicing/source/code_package/generate_all_neutral_data_add_features.r $desp $archive");
print "Post R 3 \n";
print LOG "<p>$status</p>\n";

#`cat weka.head.inner $desp.weka.data.inner.pred.body > $desp.weka.data.inner.pred.arff`;
#`cat weka.head.ss    $desp.weka.data.ss.pred.body    > $desp.weka.data.ss.pred.arff`;

#mkdir $desp if !-e $desp;
#`mv $desp.* $archive`;
print "hello\n";



my $var_c = 0;
$var_c = `wc -l < $archive/$desp.weka.data.ss.info`;
chomp($var_c);
print LOG "<p>lines of weka.data.ss.info $var_c</p>\n";
if($var_c >0){
#	print LOG "./regsnps-splicing/source/code_package/model/roc.test.sh ./regsnps-splicing/source/code_package/model/ss.4.model $archive/$desp.weka.data.ss.pred.arff $desp.splicingSite $archive 2>&1\n";
	$status =  system("~/regsnps-splicing/source/code_package/model/roc.test.sh ./regsnps-splicing/source/code_package/model/ss.4.model $archive/$desp.weka.data.ss.pred.arff $desp.splicingSite $archive 2>&1");

	print LOG "<p>$status</p>\n";
}

$var_c = `wc -l < $archive/$desp.weka.data.inner.info`;
chomp($var_c);
print LOG "<p>lines of weka.data.inner.info  $var_c</p>\n";
if($var_c >0){
	print LOG "made it here\n";
	$status = system("~/regsnps-splicing/source/code_package/model/roc.test_neutral.sh ./regsnps-splicing/source/code_package/model/combined_model_4 $archive/$desp.weka.data.inner.pred.arff $desp.exonBody $archive 2>&1");
	#$status = `./regsnps-splicing/source/code_package/model/roc.test_neutral.sh ./regsnps-splicing/source/code_package/model/exonbody.5.model $archive/$desp.weka.data.inner.pred.arff $desp.exonBody $archive 2>&1`;
#	print LOG "./regsnps-splicing/source/code_package/model/roc.test_neutral.sh ./regsnps-splicing/source/code_package/model/combined_model_4 $archive/$desp.weka.data.inner.pred.arff $desp.exonBody $archive 2>&1\n";
	print LOG "<p>$status</p>\n";

}





exit();
