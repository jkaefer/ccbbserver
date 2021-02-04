#!usr/bin/perl
use warnings;
use strict;
my $desp = $ARGV[0];
my $job = $ARGV[1];
#open(OUT, ">$job/rnafold")||die "$!";

`/N/u/jkaefer/BigRed3/.conda/envs/py37/bin/RNAfold --noPS < $job/$desp.snp.ref.mut.seq.ref.ranfold.input  > $job/$desp.snp.ref.mut.seq.ref.RNAfold`;#print OUT "$status\n";
`/N/u/jkaefer/BigRed3/.conda/envs/py37/bin/RNAfold --noPS < $job/$desp.snp.ref.mut.seq.mut.ranfold.input  > $job/$desp.snp.ref.mut.seq.mut.RNAfold`;#print OUT "$status\n";
`perl ./regsnps-splicing/source/code_package/RNAstructure/combine_data_for_RNAdistance.pl $job/$desp.snp.ref.mut.seq.ref.RNAfold $job/$desp.snp.ref.mut.seq.mut.RNAfold > $job/$desp.combined.Struct`;#print OUT "$status\n";
`/N/u/jkaefer/BigRed3/.conda/envs/py37/bin/RNAdistance < $job/$desp.combined.Struct  > $job/$desp.combined.Struct.compare`;#print OUT "$status\n";
my $cmd = "less $job/$desp.combined.Struct.compare | grep -P \'f:\\s[\\d]+\' | awk \'{print \$2}\' > $job/$desp.combined.Struct.compare.score";#print OUT "$status\n";
system($cmd);#print OUT "$status\n";
#close OUT;
