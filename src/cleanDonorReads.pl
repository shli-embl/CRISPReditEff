#usage: perl cleanDonorReads.pl input.sam output.sam chromsome donor_start donor_end
use Getopt::Long;
use strict;
use constant FALSE => 1==0;
use constant TRUE => not FALSE;
my $display_help;
my $input_sam;
my $output_sam;
my $chromosome;
my $donor_start;
my $donor_end;
my $flexible_bases=6;
## flexible_bases: interval to test is slightly extended relative to the donor, as the construct has a change to match the sequence of genome. 6 bases by defalt are considered as the maximum number of match by chance. 
my %chr_name_switch=("I"=>1,"II"=>2,"III"=>3,"IV"=>4,"V"=>5,"VI"=>6,"VII"=>7,
					"VIII"=>8,"IX"=>9,"X"=>10,"XI"=>11,"XII"=>12,"XIII"=>13,
					"XIV"=>14,"XV"=>15,"XVI"=>16);
GetOptions('h' => \$display_help, 'in=s' => \$input_sam, 'out=s' =>\$output_sam, 'chr=s' =>\$chromosome, 'd_start=s' =>\$donor_start, 'd_end=s' =>\$donor_end);
if($display_help)
{
	print "Command: \n\tperl cleanDonorReads.pl [-h] -in INPUT_SAM -out OUTPUT_SAM -chr CHROMOSOME -d_start DONOR_START_COORD -d_end DONOR_END_COORD\n\n";
	print "Function: \n\tExclude DNA donor-related reads from target sites in a given SAM file.\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput SAM.\n";
	print "\t-out\tOutput SAM (with clean target mapping).\n";
	print "\t-chr\tChromosome of target.\n";
	print "\t-d_start\tPosition of donor start on chromosome.\n";
	print "\t-d_end\tPosition of donor end on chromosome.\n";
	exit 1;
}

if((!$input_sam)||(!$output_sam)||(!$chromosome)||(!$donor_start)||(!$donor_end))
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
#my %cassette_elements=("chr15"=>"721500-723000","chr5"=>"116000-117000");

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

sub is_transloca

print "START searching for construct reads...\n";
open(INSAM,$input_sam);
my %read_is_overlapped;
my %fragment_is_contained;
my %is_softclipped_near_junction;
my %is_translocation_to_REDI;

my %discard_list_unique;
my %on_target_list_unique;
my %translocation_list_unique;
my $line_count=0;
while(my $read_line=<INSAM>)
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
		
		if($seq_self_chr eq $chromosome)
		{
			if(&is_overlapped($donor_start,$donor_end,$seq_self_map_start,$seq_self_map_start+&length_CIGAR($seq_CIGAR))
			{
				
			}
		}
		
		
		
		
		
		my $this_seq_overlap_target=FALSE;
		my $mate_seq_not_same_location=FALSE;
		#if the read overlapps target interval, returns $this_seq_overlap_target=TRUE
		if($seq_self_chr eq $chromosome)
		{
			if(($donor_start>=$seq_self_map_start)&&($donor_start<($seq_self_map_start+&length_CIGAR($seq_CIGAR)))||(($donor_end>=$seq_self_map_start)&&($donor_end<($seq_self_map_start+&length_CIGAR($seq_CIGAR))))||(($donor_start<$seq_self_map_start)&&($donor_end>=($seq_self_map_start+&length_CIGAR($seq_CIGAR)))))
			{
				$this_seq_overlap_target=TRUE;				
			}
		}
		#if the read mate maps to another chromosome, or has a fragment size > 1000, returns mate_seq_not_same_location=TRUE
		if(($seq_mate_chr ne "=")||(abs($seq_mate_map_start-$seq_self_map_start)>1000)||($fragment_size==0))
		{
			$mate_seq_not_same_location=TRUE;
		}
		#if the read overlaps target interval and its mate pair doesn't map to same loc; then check if it overlaps the junction of donor (see if it is more likely a donor read) (allows 3bp outside being aligned)
		if($this_seq_overlap_target)
		{
			if($mate_seq_not_same_location)
			{
				if((($donor_start-$seq_self_map_start)<=5)&&(($donor_end-$seq_self_map_start-&length_CIGAR($seq_CIGAR))>=(-5)))
				{
					$discard_list_unique{$seqpair_name}=1;
					
				}
				else
				{					
					#need update: if mate read map to HIS3, $is_cassette_read=1
					
					$on_target_list_unique{$seqpair_name}=1;
					
					my $mate_chr=$seq_self_chr;
					if($seq_mate_chr ne "=")
					{
						$mate_chr=$seq_mate_chr;
					}
					my $mate_coord=$seq_mate_map_start;
					
					if($cassette_elements{$mate_chr})
					{
						my @element_interval=split("-",$cassette_elements{$mate_chr});
						if(($mate_coord >= $element_interval[0])&&($mate_coord <= $element_interval[1]))
						{
							$translocation_list_unique{$seqpair_name}=1;
						}
					}
					
				}
			}
			#else if normal read pairs, check if the fragment size encompass the junction of donor
			else
			{
				my $fragment_start=&min($seq_self_map_start,$seq_mate_map_start);
				my $fragment_end=$fragment_start+$fragment_size-1;
				if((($donor_start-$fragment_start)<=6)&&($donor_end-$fragment_end>=(-6)))
				{
					$discard_list_unique{$seqpair_name}=1;
				}
				else
				{
					$on_target_list_unique{$seqpair_name}=1;
					if(&is_from_cassette($seq_CIGAR,$seq_self_map_start,$donor_start,$donor_end))
					{
						$translocation_list_unique{$seqpair_name}=1;
					}
				}
			}
		}

		
		
	}
}
close INSAM;
print "Start filtering SAM file\n";
open(INSAM,$input_sam);
open(OUTSAM,">".$output_sam);
open(DISCARD,">".$output_sam.".discarded.reads");
open(MAP,">".$output_sam.".mapped.reads");
open(TRANSLOCATION,">".$output_sam.".translocation.reads");
$line_count=0;
while(my $read_line=<INSAM>)
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
		my $FLAG=$read_column[1];
		if($discard_list_unique{$readpair_name})
		{
			print DISCARD $read_line,"\n";
		}
		else
		{
			print OUTSAM $read_line,"\n";
			if($on_target_list_unique{$readpair_name})
			{
				print MAP $read_line,"\n";
				if($translocation_list_unique{$readpair_name})
				{
					print TRANSLOCATION $read_line,"\n";
				}
			}
		}			
	}
	else
	{
		print OUTSAM $read_line,"\n";
	}
}

my $length_filtered=keys %discard_list_unique;
my $length_remain=keys %on_target_list_unique;

print "In total, ",$length_filtered," pairs are potentially from construct; ",$length_remain," pairs are overlapping the target site.\n"
