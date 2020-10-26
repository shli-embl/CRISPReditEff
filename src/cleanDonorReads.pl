#usage: perl cleanDonorReads.pl input.sam output.sam chromsome donor_start donor_end
use Getopt::Long;
use File::Basename;
use strict;
use constant FALSE => 1==0;
use constant TRUE => not FALSE;
my $display_help;
my $input_bam;
my $output_bam;
my $chromosome;
my $donor_start;
my $donor_end;
my $flexible_bases=6;
my $verbose="";
my $bl_region_file=dirname("$0")."//../data/REDI_reads_mapped_regions.txt";
## flexible_bases: interval to test is slightly extended relative to the donor, as the construct has a change to match the sequence of genome. 6 bases by defalt are considered as the maximum number of match by chance. 
my %chr_name_switch=("I"=>1,"II"=>2,"III"=>3,"IV"=>4,"V"=>5,"VI"=>6,"VII"=>7,
					"VIII"=>8,"IX"=>9,"X"=>10,"XI"=>11,"XII"=>12,"XIII"=>13,
					"XIV"=>14,"XV"=>15,"XVI"=>16);
GetOptions('h' => \$display_help, 'in=s' => \$input_bam, 'out=s' =>\$output_bam, 'chr=s' =>\$chromosome, 'd_start=i' =>\$donor_start, 'd_end=i' =>\$donor_end, 'verbose' =>\$verbose);
if($display_help)
{
	print "Command: \n\tperl cleanDonorReads.pl [-h] -in INPUT_BAM -out OUTPUT_BAM -chr CHROMOSOME -d_start DONOR_START_COORD -d_end DONOR_END_COORD\n\n";
	print "Function: \n\tExclude DNA donor-related reads from target sites in a given BAM file.\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput BAM.\n";
	print "\t-out\tOutput BAM (with clean target mapping).\n";
	print "\t-chr\tChromosome of target.\n";
	print "\t-d_start\tPosition of donor start on chromosome.\n";
	print "\t-d_end\tPosition of donor end on chromosome.\n";
	print "\t-verbose\tTurn on verbose mode.\n";
	exit 1;
}

if((!$input_bam)||(!$output_bam)||(!$chromosome)||(!$donor_start)||(!$donor_end))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

my $chr_test=$chromosome;
$chr_test=~s/chr//;
if($chr_name_switch{$chr_test})
{
	$chromosome="chr".$chr_name_switch{$chr_test};
}

sub length_CIGAR
{
	my $CIGARstring=$_[0];
	my $length_CIGAR=0;
	while($CIGARstring)
	{
		if($CIGARstring=~/^(\d+)([IDMSH])/)
		{
			my $tmp_length=$1;
			my $tmp_type=$2;
			if(($tmp_type eq "M")||($tmp_type eq "D"))
			{
				$length_CIGAR+=$tmp_length;
			}
			$CIGARstring=~s/^\d+[IDMSH]//;
		}
		elsif($CIGARstring eq "*")
		{
			$CIGARstring=~s/\*//;
		}
		else
		{
			print "CIGAR string unrecognized:",$CIGARstring,"\n";
			exit;
		}
	}
	return $length_CIGAR;
}

sub min
{
	my ($num1,$num2)=@_;
	if($num1 < $num2)
	{
		return $num1;
	}
	else
	{
		return $num2;
	}
}

sub is_overlapped
{
	my ($interval1_left,$interval1_right,$interval2_left,$interval2_right)=@_;
	if (($interval1_right < $interval2_left)||($interval1_left > $interval2_right))
	{
		return FALSE;
	}
	else
	{
		return TRUE;
	}
}

sub is_contained
#test if interval B fully contains inteval A
{
	my ($interval1_left,$interval1_right,$interval2_left,$interval2_right)=@_;
	if(($interval2_left<=$interval1_left)&&($interval2_right>=$interval1_right))
	{
		return TRUE;
	}
	else
	{
		return FALSE;
	}
}

sub is_softclipped_near_junction
#test if a read has soft clip postion matching donor boarder
#for sorting out potential REDI-target translocating reads
{
	my ($interval_left,$interval_right, $read_map_start, $CIGAR, $flexible_boundary) = @_;
	#flexible_boundary: similar as flexible bases. A softclip occurs N bases close to the donor boundary is considered as near junction, by default $flexible_boundary = 6
	my $ISNJ=FALSE;
	if($CIGAR=~/^\d+S/)
	{
		if(abs($read_map_start-$interval_left)<=$flexible_boundary)
		{
			$ISNJ=TRUE;
		}
	}
	elsif($CIGAR=~/\d+S$/)
	{
		if(abs($read_map_start+&length_CIGAR($CIGAR)-$interval_right)<=6)
		{
			$ISNJ=TRUE;
		}
	}
	return $ISNJ;
}

sub is_mapped_to_cassete
{
	my ($chr,$coord,$black_list_ref)=@_;
	my $is_mapped_to_cassete=FALSE;
	foreach my $blchr(keys %{$black_list_ref})
	{
		if($blchr eq $chr)
		{
			foreach my $blstart(keys %{$black_list_ref->{$blchr}})
			{
				foreach my $blend(keys %{$black_list_ref->{$blchr}{$blstart}})
				{
					if(($coord < $blend)&&($coord > $blstart))
					{
						$is_mapped_to_cassete=TRUE
					}
				}
			}
		}
	}
	return $is_mapped_to_cassete;
}

#Load "black-listed" regions where REDI cassete reads map to
open(BL,$bl_region_file);
my $REDI_regions;
while(my $read_line=<BL>)
{
	chomp $read_line;
	my @read_column=split("\t",$read_line);
	$REDI_regions->{$read_column[0]}{$read_column[1]}{$read_column[2]}++;
}


print "START marking donor reads...\n";
open(INBAM,"samtools view -h ".$input_bam." |");

my %overlapped_fragments;
my %discarded_fragments;
my %included_fragments;
my %REDI_translocation_fragments;
my $line_count=0;

#stopped here

while(my $read_line=<INBAM>)
{
	chomp $read_line;
	if($line_count%100000==1)
	{
		print "Processed ",$line_count-1," lines.\n";		
	}
	$line_count++;
	if($read_line!~/^\@/)
	{
		my @read_column=split("\t",$read_line);
		my ($seqpair_name,$FLAG,$seq_self_chr,$seq_self_map_start,$variable4,$seq_CIGAR,$seq_mate_chr,$seq_mate_map_start,$fragment_size)=@read_column;
		
		if(($seq_self_chr eq $chromosome)&&(&is_overlapped($donor_start,$donor_end,$seq_self_map_start,$seq_self_map_start+&length_CIGAR($seq_CIGAR))))
		#if a read overlaps target region
		{
			$overlapped_fragments{$seqpair_name}++;
			if(($seq_mate_chr eq "=")&&(abs($fragment_size)<1000))
			#if its mate read map to a near location
			{
				my $fragment_start=&min($seq_self_map_start,$seq_mate_map_start);
				my $fragment_end=$fragment_start + abs($fragment_size) - 1;
				
				if(&is_contained($fragment_start,$fragment_end,$donor_start-$flexible_bases,$donor_end+$flexible_bases))
				#judge if fragment is contained by donor boundary. if so, mark and discard the pair from target site. 
				{
					$discarded_fragments{$seqpair_name}++;
				}
				#otherwise, test if a softclip is detect near donor boundary: yes->keep the reads and mark as potential REDI translocation; no->keep the reads
				else
				{
					$included_fragments{$seqpair_name}++;
					if(&is_softclipped_near_junction($donor_start,$donor_end,$seq_self_map_start,$seq_CIGAR,$flexible_bases))
					{
						$REDI_translocation_fragments{$seqpair_name}++;
					}
				}
			}
			elsif($seq_mate_chr ne "*")
			#elseif its mate read map to a distal location(or other chromosomes)
			{
				if($seq_mate_chr eq "=")
				{
					$seq_mate_chr = $seq_self_chr;
				}
				if(&is_mapped_to_cassete($seq_mate_chr,$seq_mate_map_start,$REDI_regions))
				#test if the mate read map to REDI-related regions (from loaded black list), if so, marked and discard the pair from target site. 
				{
					if(!&is_contained($seq_self_map_start,$seq_self_map_start+&length_CIGAR($seq_CIGAR)-1,$donor_start-$flexible_bases,$donor_end+$flexible_bases))
					#if a read exceed the donor boundary (is not contained) and its mate read map to REDI-related regions, this suggests a potential REDI-target translocation. The read is kept and recorded as "potential REDI translocation"
					{
						$included_fragments{$seqpair_name}++;
						$REDI_translocation_fragments{$seqpair_name}++;
					}
					else
					#otherwise it suggests an orgin from REDI cassette. Reads are removed. 
					{
						$discarded_fragments{$seqpair_name}++;
					}
				}
				else
				#if not, keep the translocation reads. 
				{
					$included_fragments{$seqpair_name}++;
				}
			}
			else
			#elseif its mate read isn't properly mapped, discard the read
			{
				$discarded_fragments{$seqpair_name}++;
			}
		}
		
	}
}		

close INBAM;
print "Start filtering BAM file\n";
open(INBAM,"samtools view -h ".$input_bam." |");
open(OUTSAM,">".$output_bam.".tmp.sam");
open(ANNOTATION,">".$output_bam.".target_filtering.anno");
$line_count=0;

print ANNOTATION "overlapped_frag\tdiscarded_frag\tincluded_frag\tREDI_translocation_frag\n";
#output numbers of 1) all overlapped fragments to target interval; 2) excluded fragments (barcode locus-related); 3) included correctly mapped target reads; 4) potential REDI-barcode locus translocation fragments;
print ANNOTATION (my $size1 = keys %overlapped_fragments),"\t",(my $size2=keys %discarded_fragments),"\t",(my $size3=keys %included_fragments),"\t",(my $size4=keys %REDI_translocation_fragments),"\n";
if($verbose)
{
	print ANNOTATION "\n\nREAD_INFO:\n";
}

while(my $read_line=<INBAM>)
{
	chomp $read_line;
	if($line_count%100000==1)
	{
		print "Processed ",$line_count-1," lines.\n";
	}
	$line_count++;
	if($read_line!~/^\@/)
	{
		my @read_column=split("\t",$read_line);
		my $readpair_name=$read_column[0];
		if($overlapped_fragments{$readpair_name})
		{
			if($discarded_fragments{$readpair_name})
			{
				if($verbose)
				{
					print ANNOTATION "DISCARD: ",$readpair_name,"\n";
				}
				
			}
			elsif($included_fragments{$readpair_name})
			{
				print OUTSAM $read_line,"\n";
				print ANNOTATION "INCLUDE: ",$readpair_name,"\n";
				if($REDI_translocation_fragments{$readpair_name})
				{
					if($verbose)
					{
						print ANNOTATION "REDITRANS: ",$readpair_name,"\n";
					}
					
				}
			}
			else
			{
				
				if($verbose)
				{
					print ANNOTATION "UNKOWN: ",$readpair_name,"\n";
				}
				print $read_line,"\n";
				
			}
		}
		else
		{
			print OUTSAM $read_line,"\n";
		}
	}
	else
	{
		print OUTSAM $read_line,"\n";
	}
}

system("samtools view -bS ".$output_bam.".tmp.sam > ".$output_bam);
system("samtools index $output_bam");
system("rm ".$output_bam.".tmp.sam");
