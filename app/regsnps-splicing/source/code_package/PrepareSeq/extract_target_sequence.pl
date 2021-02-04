#!usr/bin/perl
use warnings;

## generate the ref sequence and mutated sequence for each neutral SNP
my $snpfile = "";
my $desp = "";
my $job = "";
if ($#ARGV == 2){
	$snpfile = $ARGV[0];
	$desp    = $ARGV[1];
	$job     = $ARGV[2];
}else{
	print "Usage:perl PrepareSeq/extract_target_sequence.pl Input short_description_of_input\n\n\n";
	exit();
}


#my $snpfile = 'test';
my $mcfile = "/N/slate/jkaefer/splicingdb/hg19.exon.300bp.extension.no-alt-anno.unique";
#my $maf_snp = 'snp_maf';
#open(SM, ">$maf_snp")||die "$!";
open(LOG, ">$job/log");
open (IN, $snpfile)||die "$!";
open(CP, ">$snpfile.cp")||die "$!";
while(<IN>){
	s/\s+$//;
	s/^chr//;
	print CP "$_\n";
}
close CP;
close IN;
`mv $snpfile.cp $snpfile`;

my %mut = ();
open (IN, $snpfile)||die "$!";
while(<IN>){
	#print LOG "Hello\n";
	s/\s+$//;
	if(!/^#/){
		s/^chr//;
		my @line = split(/\t|\s+/, $_);
		my $chr = "";
		my $pos = "";
		my $id = "";
		my $refa = "";
		my @alta = "";
		my $maf =  0;

		if($_ =~ /^([0-9XYxy]+):(\d+)\s+([ACGT])\s+([ACGT])$/){#----------- general format: chr22:18020382	ref	mut
		#if ($#line eq 2){#----------- general format: chr22:18020382    ref     mut
			print LOG "gen form\n";
			$chr = $1;
			$pos = $2;
			$id = $line[0];
			$refa = $line[1];
			@alta = split(',', $line[2]);
			$maf = 1;
		}elsif (/^(\d+)\s+(\d+)\s+(\d+)\s+([ACGT])\s+([ACGT])/){# VCF format
			print LOG "vcf form\n";
			$chr = $line[0];
			$pos = $line[1];
			$id = "chr$line[0]:$line[1]";
			$refa = $line[3];
			@alta = split(',', $line[4]);
			$line[$#line] =~ /;AF=([0-9\.]+);{0,1}/;
			$maf = $1;
		}elsif(/^([0-9xyXY]+)\s+(\d+)\s+([0-9a-zA-Z\.]+)\s+([ACGT])\s+([ACGT,]+)/){## VCF 4.0, 4.1
			print LOG "new vcf\n";
			$chr = $line[0];
			$pos = $line[1];
			$id = "chr$line[0]:$line[1]:";
			$refa = $line[3];
			@alta = split(',', $line[4]);
			$line[$#line] =~ /;CAF=([0-9\.]+),([0-9\.]+);{0,1}/;
			$maf = $2;
		#}
		#elsif(/^chr([0-9xyXY]+)\s+(\d+)\s+([0-9a-zA-Z\.]+)\s+([ACGT])\s+([ACGT,]+)/){ ## vcf4.1
		#	#s/^chr//;
		#	#$_ =~ /^([0-9xyXY]+)\s+(\d+)\s+([0-9a-zA-Z\.]+)\s+([ACGT])\s+([ACGT,]+)/
		#	$line[0] =~ s/chr//;
		#	$chr = $line[0];
		#	#$chr =~ s/chr//;
		#	$pos = $line[1];
		#	$id = "chr$chr:$line[1]:";
		#	$refa = $line[3];
		#	@alta = split(',', $line[4]);
		#	$line[$#line] =~ /;CAF=([0-9\.]+),([0-9\.]+);{0,1}/;
		#	$maf = $2;
		}else{
			#print "Please input your data in the following format:\n1.chr22:18020382 ref mut\n2.VCF format\n";
			#exit();
			print "***** NOT required format ---> \t$_\n";
			next;
		}
		foreach my $al (@alta){
			$mut{"chr$chr"}{length($pos)}{$pos}{"$refa-$al"}{'ref'} = $refa;
			$mut{"chr$chr"}{length($pos)}{$pos}{"$refa-$al"}{'mut'} = $al;
			$mut{"chr$chr"}{length($pos)}{$pos}{"$refa-$al"}{'maf'} = $maf;
		}
	}
}
close IN;

open(OUT, ">$job/$desp.snp.ref.mut.seq")||die "$!";
open(PRO, ">$job/$desp.snp.ref.mut.seq.proximity")||die "$!";
open(EVO, ">$job/$desp.snp.ref.mut.seq.mutBed")||die "$!"; # -5 to +5 pos around the SNP
open(WRN, ">$job/$desp.invalid_input")||die "$!";




my %revc = ();
$revc{'A'} = 'T';
$revc{'T'} = 'A';
$revc{'G'} = 'C';
$revc{'C'} = 'G';
$revc{'a'} = 'T';
$revc{'t'} = 'A';
$revc{'g'} = 'C';
$revc{'c'} = 'G';

open (IN, $mcfile)||die "$!";
while(<IN>)
{
	#print LOG "two\n";
	s/\s+$//;
	my @line = split(/\t/, $_);
	my @tmp = split(/:/, $line[0]);
	my @pos = ();
	$pos[0] = $tmp[0];#chr
	my @tt  = split(/-/, $tmp[1]);
	$pos[1] = $tt[0];#start
	$pos[2] = $tt[1];#end
	$pos[3] = $tmp[2];#strand
	#print LOG "position\n";
	$line[2] =~ /([atgcn]+)([ATGCN]+)([atgcn]+)/;
	my $upintron = $1;
	my $exon = $2;
	my $dnintron = $3;	
	
	my $hit = 0;
	#chromosomes at first 
	foreach my $m  (keys %{$mut{$pos[0]}{length($pos[1])}}){
		#print LOG "outer\n";
        foreach my $mt (keys %{$mut{$pos[0]}{length($pos[1])}{$m}} ){
        	#print LOG "inner\n";
		if ($m >= ($pos[1]+length($upintron)) and $m <= ($pos[2]-length($dnintron))){
			print LOG "innerinner\n";
			my $refs = $line[2];
			my $muts = $refs;
			my $offs = $m - $pos[1];
			my $nt_in_1kg = $mut{$pos[0]}{length($pos[1])}{$m}{$mt}{'ref'};# nt from %mut is always from pos strand
			my $alt_nt    = $mut{$pos[0]}{length($pos[1])}{$m}{$mt}{'mut'};
			my $rfvar = $nt_in_1kg;
	
			my $pro5 = $offs - 300;;
			my $pro3 = $pos[2]-$m-300;
			#flipping based on strand
			if ($pos[3] eq '-'){
				print LOG "neg\n";
				$offs = $pos[2]-$m;
				$nt_in_1kg = $revc{$nt_in_1kg};
				$alt_nt    = $revc{$alt_nt};
				
				my $tmp = $pro5;
				$pro5 = $pro3;
				$pro3 = $tmp;
			}
			if ( substr($refs, $offs, 1) eq $nt_in_1kg ){ # just make sure they are the same
				print LOG "same\n";
				substr($muts, $offs, 1, $alt_nt );
				print OUT "$line[0]\t$line[1]\t$pos[0]:$m:$mt\t$refs\t$muts\n";
				print PRO "$line[0]\t$line[1]\t$pos[0]:$m:$mt\t$pro5\t$pro3\t$mut{$pos[0]}{length($pos[1])}{$m}{$mt}{'maf'}\n";
				my $lp = $m -5;
				my $rp = $m +5;
				print EVO "$pos[0]:$lp-$rp\n";
				$hit = 1;
			}
			else{
				my $yourvar = "$pos[0]:$m:[ref]$rfvar";
				my $refvar  = "$pos[0]:$m:[ref]".substr($refs, $offs, 1);
				print "Your variant $yourvar  does not match refseq  \n";
				print WRN "Your variant $yourvar  does not match refseq  \n";
			}
		}
	}}

	if (length($pos[1]) != length($pos[2]) and $hit == 0){
		foreach my $m  (keys %{$mut{$pos[0]}{length($pos[1])}}){
		foreach my $mt (keys %{$mut{$pos[0]}{length($pos[1])}{$m}} ){
			if ($m >= ($pos[1]+length($upintron)) and $m <= ($pos[2]-length($dnintron))){
				my $refs = $line[2];
				my $muts = $refs;
				my $offs = $m - $pos[1];
				my $nt_in_1kg = $mut{$pos[0]}{length($pos[1])}{$m}{$mt}{'ref'};# nt from %mu is always from pos strand
				my $alt_nt    = $mut{$pos[0]}{length($pos[1])}{$m}{$mt}{'mut'};
				if ($pos[3] eq '-'){
					$offs = $pos[2]-$m;
					$nt_in_1kg = $revc{$nt_in_1kg};
					$alt_nt    = $revc{$alt_nt};
				}
				if ( substr($refs, $offs, 1) eq $nt_in_1kg ){ # just make sure they are the same
					substr($muts, $offs, 1, $alt_nt );
					print OUT "$line[0]\t$line[1]\t$pos[0]:$m:$mt\t$refs\t$muts\n";
					$hit = 1;
				}
				else{
					my $yourvar =  "$pos[0]:$pos[1]:[ref]$rfvar";
					my $refvar = "$pos[0]:$pos[1]:[ref]".substr($refs, $offs, 1);
					print "Your variant $yourvar  does not match refseq  \n";
					#print "$line[0] not match\n";
					print WRN "Your variant $yourvar  does not match refseq  \n";
				}
			}
		}
	}}

}
close IN;
close OUT;
close PRO;
close EVO;
close WRN;


%mut = ();
