#!usr/bin/perl 
#use warnings;
use Fcntl;
## local the exon where a SNP is located
## and the position of that exon in protein


## PTM annotation is extracted from uniport
## Pfam annotation is added with clan information


#input format chr1	xxx	+, 0-based

my $filename = '/N/slate/jkaefer/splicingdb/hg19.protein.pfam.data'; # pfam file
#-----------------------------------------index whole genome pfam file-----------------------------
open(my $pfam_orig, "<", $filename) 
     or die "Can't open $filename for reading: $!";

# open the index and build it if necessary
# there's a race condition here: two copies of this
# program can notice there's no index for the file and
# try to build one.  This would be easily solved with
# locking
$indexname = "$filename.index";
sysopen(my $pfam_idx, $indexname, O_CREAT|O_RDWR)
			or die "Can't open $indexname for read/write: $!";
build_index($pfam_orig, $pfam_idx,'pfam') if -z $indexname;  # XXX: race unless lock

## ----------------------------------------index whole genome ptm file------------------------
#my $ptm_filename = 'uniprot_sprot.dat.human_PTM';#uniport
my $ptm_filename = '/N/slate/jkaefer/splicingdb/dbPTM3.human.human_PTM';
#my  $ptm_filename = "hg19.transcript.wholeGenome.PTM"; #modpred

open(my $ptm_orig, "<", $ptm_filename) 
     or die "Can't open $ptm_filename for reading: $!";

$ptm_indexname = "$ptm_filename.index";
sysopen(my $ptm_idx, $ptm_indexname, O_CREAT|O_RDWR)
			or die "Can't open $ptm_indexname for read/write: $!";
build_index($ptm_orig, $ptm_idx,'ptm') if -z $ptm_indexname;  # XXX: race unless lock


my %ptm_line_to_id = ();## mapping id to its line number
my %pfam_line_to_id = ();## mapping id to its line number


#my $infile = 'Disease_causing_splicing_samesense';
my $infile = $ARGV[0];
open (OUT, ">$infile.post-trans.mod")||die "$!";
open (MA,  ">$infile.pfam-matrix")||die "$!";
open (MO,  ">$infile.ptm-matrix.dbPTM")||die "$!";
open (DF,  ">$infile.disorder-matrix")||die "$!";

print MO "id\tgene";
my @clans = ();
open (CL, '/N/slate/jkaefer/splicingdb/clans')||die "$!";
print MA "id";
while(<CL>)
{
	s/\s+$//;
	#$clans{$_} = 1;
	push @clans, $_;
	print MA "\t$_";
}
close CL;
print MA "\n";

my $exon_loc_fl = '/N/slate/jkaefer/splicingdb/transcript_aaSeq.good'; # removed some transcripts which don't exist in ucsc browser data. Like NM_001242461
open(IN, "$exon_loc_fl")||die "$!";
my %genes = ();
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	my @bits_start = split('', $line[5]);
	my @bits_end   = split('', $line[6]);
	$genes{$line[1]}{$#bits_start+1}{$#bits_end+1}{$line[0]} = $_; #use chr and location information to group genes
}
close IN;


my %lo_conver_ptm = ();
open (IN, "/N/slate/jkaefer/splicingdb/ptm.line_to_offset_convertor")||die "$!";# gene may have multiple lines
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	@{$lo_conver_ptm{$line[0]}} = @line[1..2];
    push @{$ptm_line_to_id{$line[3]}} , $line[0];
}
close IN;


my %lo_conver_pfam = ();
open (IN, "/N/slate/jkaefer/splicingdb/pfam.line_to_offset_convertor")||die "$!";
while(<IN>)
{
	s/\s+$//;
	my @line = split(/\t/, $_);
	@{$lo_conver_pfam{$line[0]}} = @line[1..2];
    push @{$pfam_line_to_id{$line[3]}} , $line[0];
}
close IN;

my %ptm_table = ();
my @ptm_name = ();
open(IN, "/N/slate/jkaefer/splicingdb/ptm.list")||die "$!";#uniport
#open(IN, "modpred.ptm.list")||die "$!"; #modpred
while(<IN>)
{
	s/\s+$//;
	$ptm_table{$_} = 0;
	push @ptm_name, $_;
	print MO "\t$_";
}
close IN;
print MO "\n";

######################### process each input SNP ############################
open (IN, $infile) ||die "$!"; # HGMD/1000G input is 1-based, the hg19.bed is 0-based
while(<IN>)
{
	s/\s+$//;
	my @tmp = split(/\t/, $_);
	my $snp = $tmp[2];

	#my $snp = $_;
	my @pos = split(/:/, $snp);
	$pos[1] -= 1;

	my $loc_width = length($pos[1]);
	my $loc_first_num = substr($pos[1],0,1);

	my @exon_loc = ();
	#@exon_loc = snp_loc(\$genes{$pos[0]}{$loc_width}{$loc_width}, $snp );
	@exon_loc = snp_loc(\$genes{$pos[0]}{$loc_width}{$loc_width}, \@tmp ); # when in put is hgmd.exonic.snp.ref.mut.seq 
	print OUT "$snp\t$exon_loc[0]:$exon_loc[1]:$exon_loc[2]\t";
	if ($#exon_loc eq -1)
	{
		#@exon_loc = snp_loc(\$genes{$pos[0]}{$loc_width}{$loc_width+1}, $snp );
		@exon_loc = snp_loc(\$genes{$pos[0]}{$loc_width}{$loc_width}, \@tmp );
		print OUT "$exon_loc[0]:$exon_loc[1]:$exon_loc[2]\t";
		if ($#exon_loc eq -1)
		{
			#@exon_loc = snp_loc( \$genes{$pos[0]}{$loc_width-1}{$loc_width}, $snp );
			@exon_loc = snp_loc(\$genes{$pos[0]}{$loc_width}{$loc_width}, \@tmp );
			print OUT "$exon_loc[0]:$exon_loc[1]:$exon_loc[2]\t";
		}
	}
	
	if ($#exon_loc > 1)
	{
		my $line_number = 'hg19_refGene_'.$exon_loc[0];#gene id		
		my @pfam_lines = @{$pfam_line_to_id{$line_number}};
		
		$line_number =~ s/^hg19_refGene_//; #uniprot
		my @ptm_lines  = @{$ptm_line_to_id{$line_number}};# use the first record
		
		#-----------------------------[ extract Pfam information ]
		my $pfam_c = 0;
		my %clan_overlap = ();
		my @pfam_info = ();
		foreach my $ln (@pfam_lines)
		{
			my $line = line_with_index($pfam_orig, $pfam_idx, $ln+1,\%lo_conver_pfam);
			if ($line ne ''){
				#print "$snp\t$exon_loc[0]:$exon_loc[1]:$exon_loc[2]\t";
				@pfam_info = split(/\s+/, $line);
				if ($pfam_info[1] >= $exon_loc[1] and $pfam_info[1] <= $exon_loc[2] or
					$pfam_info[2] >= $exon_loc[1] and $pfam_info[2] <= $exon_loc[2] or
					$pfam_info[1] <= $exon_loc[1] and $pfam_info[2] >= $exon_loc[2] or
					$pfam_info[1] >= $exon_loc[1] and $pfam_info[2] <= $exon_loc[2])
					{
						my $overl = 0;
						my $a = $pfam_info[1];
						my $b = $pfam_info[2];
						my $c = $exon_loc[1];
						my $d = $exon_loc[2];
						if ($a >= $c and $a <= $d and $d <= $b){
							$overl = $d - $a+1;
						}elsif ($c >= $a and $c <= $b and $d >= $b){
							$overl = $b-$c+1;
						}elsif ($c >$a and $d < $b){
							$overl = $d-$c+1;
						}elsif($a>$c and $b < $d){
							$overl = $b-$a+1;
						}else{
							$overl=-1;
						}
						
						$overl /= ($b-$a+1); # modified on Dec,1,2014
					
						$overl = sprintf("%.3f", $overl);
						print OUT "$pfam_info[1] $pfam_info[2] $pfam_info[5] $pfam_info[6] $pfam_info[10] $pfam_info[14] $overl;";
						if ($pfam_info[14] eq 'No_clan'){
							$clan_overlap{$pfam_info[14]} .= "$pfam_info[6]:$overl;";
							#if ($clan_overlap{$pfam_info[14]} ne ''){
							#	if ($overl > $clan_overlap{$pfam_info[14]}){
							#		$clan_overlap{$pfam_info[14]} = $overl;
							#	}
							#}else{
							#	$clan_overlap{$pfam_info[14]} .= "$overl";
							#}
						}
						else{
							$clan_overlap{$pfam_info[14]} .= "$overl;";
							#if ($clan_overlap{$pfam_info[14]} ne ''){
							#	if ($overl > $clan_overlap{$pfam_info[14]}){
							#		$clan_overlap{$pfam_info[14]} = $overl;
							#	}
							#}else{
							#	$clan_overlap{$pfam_info[14]} = "$overl";
							#}
						}
						$pfam_c += 1;
					}
			}
			#else
			#{
			#	print OUT "no location is found for this snp\n";
			#}
		
			#print OUT "\t";
		}
		print MA "$snp";
		foreach my $c (@clans)
		{
			if ($clan_overlap{$c} eq '')
			{
				print MA "\t0";
			}
			else
			{
				print MA "\t$clan_overlap{$c}";
			}
		}
		print MA "\n";

		print OUT "|";
		#--------------[ extracted PTM information ]
		my $ptm_c = 0;
		
		if ($#ptm_lines >=0)
		{
			my $line = line_with_index($ptm_orig, $ptm_idx, $ptm_lines[0]+1, \%lo_conver_ptm);
			#die "Didn't find line $line_number in $filename" unless defined $line;
			#print $line;
			if ($line ne '')
			{
				my @ptm_sites = split(/[\t\|]/, $line); #uniport
				#my @ptm_sites = split(/[\t;]/, $line);   #modpred
				for my $k (2..($#ptm_sites-1)) #uniport
				#for my $k (1..($#ptm_sites-1))  #modpred
				{
					#$ptm_sites[$k]=~ /^[A-Z]{1}([0-9]+)/;      # modpred
					my @ptm_info = split(/;/, $ptm_sites[$k]);# uniport
					#my @ptm_info = split(/\s/, $ptm_sites[$k]);# modpred
					my $aa_pos = $ptm_info[0]; #uniport
					#my $aa_pos = $1;            #modpred

					if ($aa_pos >= $exon_loc[1] and $aa_pos <= $exon_loc[2])
					{
						print OUT "$ptm_sites[$k]\@";
						$ptm_c +=1;
						$ptm_table{$ptm_info[1]} +=1;
					}
				}
			}
			#print OUT "\t";
			$ptm_c /= ($exon_loc[2] - $exon_loc[1]+1);
		}
		print OUT "|";
		print MO "$snp\t$exon_loc[0]:$exon_loc[1]:$exon_loc[2]";
		foreach my $p (@ptm_name)
		{
			#normalize PTM count ## Dec,1,2014
			#print MO "\t$ptm_table{$p}";# comment on Dec,1,2014
			my $norm_ptm = $ptm_table{$p}/(($exon_loc[2] - $exon_loc[1]+1)/300); # modified on Dec,1,2014
			print MO "\t$norm_ptm";
			$ptm_table{$p} =0;
		}
		print MO "\n";

		#----------------------[ extract disorder score information]
		my $dsfile  = "/N/slate/jkaefer/splicingdb/disorder_score/hg19_refGene_".$exon_loc[0].".meta";
		my $ave_dis = 0;
		my @disorder_sites = ();
		my $dis_feature = "";
		if (-f $dsfile )
		{
			my $lc = 0;
			open (DS, "$dsfile")||die "$!";
			while(<DS>)
			{
				s/\s+$//;
				s/^\s+//;
				$lc +=1;
				my @ds = split(/\s+/, $_);
				if ($lc >= $exon_loc[1] and $lc <= $exon_loc[2])
				{
					print OUT "$ds[1] $ds[2];";
					$ave_dis += $ds[2];
					push @disorder_sites, $ds[2];

				}
			}
			close DS;
			$dis_feature = disorder_structure(\@disorder_sites);
			$ave_dis /= ($exon_loc[2] - $exon_loc[1]  +1);
			#print OUT "\n";
		}
		print OUT "|\t";
		print OUT "$pfam_c\t$ptm_c\t$ave_dis\n";
		if($dis_feature ne ''){
			print DF "$snp\t$exon_loc[0]:$exon_loc[1]:$exon_loc[2]\t$dis_feature\n";
		}else{
			#print "";
		}
	}
	else
	{
		print OUT "\n";
		print "no location is found for this snp $snp\t Maybe located in UTR or introns\n";
	}
}
close IN;
close MA;
close OUT;
close MO;

sub disorder_structure
{
	my $ds = shift;

	my $min = 1;
	my $max = 0;
	my $ave_dis = 0;
	my $ave_ord = 0;
	my $ave_score = 0;
	my $dis_c = 0;
	my @order_segments = (); #length of each segment
	my @disor_segments = (); #length of each segment
	my $switch = 0; # number of switches between ordered/disordered regions
	my $order_len = 0;
	my $disor_len = 0;

	foreach my $d ( 0..$#{$ds} ){
		$ave_score += $$ds[$d];
		if ($$ds[$d] > $max){
			$max = $$ds[$d];
		}
		if ($$ds[$d] < $min){
			$min = $$ds[$d];
		}
		if ($$ds[$d] > 0.5)	{
			$ave_dis += $$ds[$d];
			$dis_c +=1;
		}else{
			$ave_ord += $$ds[$d];
		}
		
		if ($d >0){
			if ($$ds[$d] > 0.5 and $$ds[$d-1] <= 0.5){
				$switch +=1;
				$order_len = sprintf("%.3f", $order_len/($#{$ds}+1));
				push @order_segments, $order_len;
				$order_len = 0;
				$disor_len = 1;
			}elsif( $$ds[$d] <= 0.5 and $$ds[$d-1] > 0.5  ){
				$switch +=1;
				$disor_len = sprintf("%.3f",$disor_len/($#{$ds}+1));
				push @disor_segments, $disor_len;
				$order_len = 1;
				$disor_len = 0;
			}
			elsif($$ds[$d] > 0.5 and $$ds[$d-1] > 0.5 ){
				$disor_len += 1;
			}
			elsif($$ds[$d] <= 0.5 and $$ds[$d-1] <= 0.5){
				$order_len += 1;
			}
		}else{
			if ($$ds[$d] > 0.5){
				$disor_len +=1;
			}else{
				$order_len +=1;
			}
		}
	}
###################################
	if ($switch == 0){
		if ($disor_len == 0){
			$order_len = sprintf("%.3f", $order_len/($#{$ds}+1));
			push @order_segments, $order_len;
			push @disor_segments, 0;
		}else{
			$disor_len = sprintf("%.3f",$disor_len/($#{$ds}+1));
			push @disor_segments, $disor_len;
			push @order_segments, 0;
		}
	}else{
		if ($$ds[$#{$ds}] > 0.5 ){
			$disor_len = sprintf("%.3f",$disor_len/($#{$ds}+1));
			push @disor_segments, $disor_len;
		}else{
			$order_len = sprintf("%.3f", $order_len/($#{$ds}+1));
			push @order_segments, $order_len;
		}
	}

	if ($dis_c != 0 and $dis_c != ($#{$ds} +1))	{
		$ave_dis /= $dis_c;
		$ave_ord /= ($#{$ds} +1 - $dis_c);
	}elsif($dis_c == 0)	{
		$ave_dis = 0;
		$ave_ord /= ($#{$ds} +1 - $dis_c);
	}elsif ($dis_c == ($#{$ds} +1))	{
		$ave_dis /= $dis_c;
		$ave_ord = 0;
	}
	$ave_dis = sprintf("%.3f", $ave_dis);
	$ave_ord = sprintf("%.3f", $ave_ord);
	
	my $dis_len_ave = 0;
	my $ord_len_ave = 0;
	@order_segments = sort {$a <=> $b} @order_segments;
	@disor_segments = sort {$a <=> $b} @disor_segments;
	foreach my $s (@order_segments){
		$ord_len_ave += $s;
	}
	foreach my $s (@disor_segments){
		$dis_len_ave += $s;
	}
	if ($#order_segments >= 0){
		$ord_len_ave /= ($#order_segments+1);
	}
	if ($#disor_segments >=0){
		$dis_len_ave /= ($#disor_segments+1);
	}

	$ave_score = $ave_score/($#{$ds}+1);
	$ord_len_ave = sprintf("%.3f", $ord_len_ave);
	$dis_len_ave = sprintf("%.3f", $dis_len_ave);
	$ave_score   = sprintf("%.3f", $ave_score);
	my $feature= "$min\t$max\t$ave_dis\t$ave_ord\t$switch\t$dis_len_ave\t$ord_len_ave\t$disor_segments[0]\t$disor_segments[-1]\t$order_segments[0]\t$order_segments[-1]\t$ave_score";
	return ($feature);
	#@order_segments @disor_segments

}


sub snp_loc
{
	my $data = shift;
	my $snp = shift;
	
	my @pos = split(/:/, $$snp[2]);
	$pos[1] -=1;
	my $gname = $$snp[1];
	$gname =~ s/_\d+$//;

	foreach my $g (keys %$$data){
		my @ge = split(/\t/, $$$data{$g});
		#if ($pos[1] >= $ge[5] and $pos[1] <= $ge[6] and $pos[2] eq $ge[2] and $g eq $gname ) # the same strand
		if ($pos[1] >= $ge[5] and $pos[1] <= $ge[6]  and $g eq $gname )
		{
			my @exon_l = split(/,/, $ge[8]);
			my @exon_r = split(/,/, $ge[9]);
			my @exon_len = split(/,/, $ge[12]);
			for my $i (0..$#exon_l){
				if ($pos[1] >= $exon_l[$i] and $pos[1] <= $exon_r[$i]){
					my $start = 0;
					my $end   = 0;
					
					if ($ge[2] eq '+'){
						for my $j (0..$i){
							$start += $exon_len[$j];
							$end   += $exon_len[$j];
						}
						$start -= $exon_len[$i]; # nucleotide location
					}
					elsif($ge[2] eq '-'){
						for my $j (0..($#exon_len-$i)){
							$start += $exon_len[$j];
							$end   += $exon_len[$j];
						}
						$start -= $exon_len[$#exon_len-$i]; # nucleotide location
					}
					else{}
					
					if($start % 3){
						$start = int($start/3);
					}else{
						$start = int($start/3) +1; # first amino acid location is 1
					}
					if($end % 3){
						$end   = int($end/3);
					}else{
						$end   = int($end/3)   +1;
					}
					
					
					#return("$ge[0]:$start:$end")
					return($ge[0],$start,$end);
				}
			}
		}
	}
}


#----------------------------------------------------------------------------------------------------
  # usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
  # returns line or undef if LINE_NUMBER was out of range
  sub line_with_index
  {
      my $data_file   = shift;
      my $index_file  = shift;
      my $line_number = shift;
      my $con         = shift;
      
      my $size;               # size of an index entry
      my $i_offset;           # offset into the index of the entry
      my $entry;              # index entry
      my $d_offset;           # offset into the data file
      
      #$size = length(pack("N", 0));
      #$i_offset = $size * ($line_number-1);
      #seek($index_file, $i_offset, 0) or return; # set handler to the wanted line
      #read($index_file, $entry, $size); # read a line
      #$d_offset = unpack("N", $entry); # unpack read content
      $d_offset = $$con{$line_number -1}[0];
      seek($data_file, $d_offset, 0);
      return scalar(<$data_file>);
  }
  
# usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
sub build_index
{
    my $data_file  = shift;
    my $index_file = shift;
    my $data_type = shift;
    my $offset     = 0;
    
    open( OF, ">$data_type.line_to_offset_convertor")||die "$!";
    my $linc = 0;
    my $len = 0;
    while (<$data_file>)
    {
        my @line = split(/\t|\s+/,$_);
        print OF "$linc\t$offset\t";
        $len = $offset;
        print $index_file pack("N", $offset);
        $offset = tell($data_file);
        $len = $offset - $len;
        print OF "$len\t$line[0]\n";
        #$$l2i{$line[0]} = $linc;
        $linc ++;
    }
    close OF;
}












