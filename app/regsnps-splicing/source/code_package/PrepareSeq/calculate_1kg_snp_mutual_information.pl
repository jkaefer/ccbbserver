#!usr/bin/perl
use warnings;

## generate the ref sequence and mutated sequence for each neutral SNP
#my $snpfile = "/data/pro/as/data/1000Genomes/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.clean.annovar.snp.exonic_variant_function";
my $snpfile = $ARGV[0];
#my $snpfile = 'test';
my $mcfile = "/N/slate/jkaefer/splicingdb/hg19.exon.300bp.extension.no-alt-anno.unique";
my $maf_snp = 'snp_maf';
open(SM, ">$maf_snp")||die "$!";

my %mut = ();
open (IN, $snpfile)||die "$!";
open (OUT,">neutral.snp.ref.mut.seq")||die "$!";
open(PRO, ">neutral.snp.ref.mut.seq.proximity")||die "$!";
open(EVO, ">neutral.snp.ref.mut.seq.mutBed")||die "$!"; # -5 to +5 pos around the SNP
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	if ($line[1] =~ /^synonymous/ or $line[1] =~ /nonsynonymous/)
#	if ($line[1] =~ /^synonymous/)
	{
		## MAF > 0.1
		$line[$#line] =~ /;AF=([0-9\.]+);{0,1}/;
		my $maf = $1;
		if ($maf >= 0.02)
		{
			my $id = "chr$line[3]:$line[4]";
			my $refa = $line[6];
			my $alta = $line[7];
			$mut{"chr$line[3]"}{length($line[4])}{$line[4]}{'ref'} = $line[6];
			$mut{"chr$line[3]"}{length($line[4])}{$line[4]}{'mut'} = $line[7];
			$mut{"chr$line[3]"}{length($line[4])}{$line[4]}{'maf'} = $maf;

			print SM "$id\t$maf\n";
		}
	}
}
close IN;
close SM;

my %revc = ();
$revc{'A'} = 'T';
$revc{'T'} = 'A';
$revc{'G'} = 'C';
$revc{'C'} = 'G';

open (IN, $mcfile)||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	my @tmp = split(/:/, $line[0]);
	my @pos = ();
	$pos[0] = $tmp[0];#chr
	my @tt  = split(/-/, $tmp[1]);
	$pos[1] = $tt[0];#start
	$pos[2] = $tt[1];#end
	$pos[3] = $tmp[2];#strand
	
	$line[2] =~ /([atgcn]+)([ATGCN]+)([atgcn]+)/;
	my $upintron = $1;
	my $exon = $2;
	my $dnintron = $3;	
	
	my $hit = 0;
	foreach my $m (keys %{$mut{$pos[0]}{length($pos[1])}})
	{
		if ($m >= ($pos[1]+length($upintron)) and $m <= ($pos[2]-length($dnintron)))
		{
			my $refs = $line[2];
			my $muts = $refs;
			my $offs = $m - $pos[1];
			my $nt_in_1kg = $mut{$pos[0]}{length($pos[1])}{$m}{'ref'};# nt from %mut is always from pos strand
			my $alt_nt    = $mut{$pos[0]}{length($pos[1])}{$m}{'mut'};
			
			my $pro5 = $offs - 300;;
			my $pro3 = $pos[2]-$m-300;

			if ($pos[3] eq '-')
			{
				$offs = $pos[2]-$m;
				$nt_in_1kg = $revc{$nt_in_1kg};
				$alt_nt    = $revc{$alt_nt};
				
				my $tmp = $pro5;
				$pro5 = $pro3;
				$pro3 = $tmp;
			}
			if ( substr($refs, $offs, 1) eq $nt_in_1kg ) # just make sure they are the same
			{
				substr($muts, $offs, 1, $alt_nt );
				print OUT "$line[0]\t$line[1]\t$pos[0]:$m\t$refs\t$muts\n";
				print PRO "$line[0]\t$line[1]\t$pos[0]:$m\t$pro5\t$pro3\t$mut{$pos[0]}{length($pos[1])}{$m}{'maf'}\n";
				my $lp = $m -5;
				my $rp = $m +5;
				print EVO "$pos[0]:$lp-$rp\n";
				$hit = 1;
			}
			else
			{
				print "$line[0] not match\n";
			}
		}
	}

	if (length($pos[1]) != length($pos[2]) and $hit == 0)
	{
		foreach my $m (keys %{$mut{$pos[0]}{length($pos[1])}})
		{
			if ($m >= ($pos[1]+length($upintron)) and $m <= ($pos[2]-length($dnintron)))
			{
				my $refs = $line[2];
				my $muts = $refs;
				my $offs = $m - $pos[1];
				my $nt_in_1kg = $mut{$pos[0]}{length($pos[1])}{$m}{'ref'};# nt from %mu is always from pos strand
				my $alt_nt    = $mut{$pos[0]}{length($pos[1])}{$m}{'mut'};
				if ($pos[3] eq '-')
				{
					$offs = $pos[2]-$m;
					$nt_in_1kg = $revc{$nt_in_1kg};
					$alt_nt    = $revc{$alt_nt};
				}
				if ( substr($refs, $offs, 1) eq $nt_in_1kg ) # just make sure they are the same
				{
					substr($muts, $offs, 1, $alt_nt );
					print OUT "$line[0]\t$line[1]\t$pos[0]:$m\t$refs\t$muts\n";
					$hit = 1;
				}
				else
				{
					print "$line[0] not match\n";
				}
			}
		}
	}

}
close IN;
close OUT;
close PRO;
close EVO;
















