#usage: perl step6_cleanDonorReads.pl input.sam output.sam chromsome donor_start donor_end
use strict;
use constant FALSE => 1==0;
use constant TRUE => not FALSE;
my $input_sam=$ARGV[0];
my $output_sam=$ARGV[1];
my $chromosome=$ARGV[2];
my $chr_test=$chromosome;
$chr_test=~s/chr//;
my $donor_start=$ARGV[3];
my $donor_end=$ARGV[4];
my %chr_name_switch=("I"=>1,"II"=>2,"III"=>3,"IV"=>4,"V"=>5,"VI"=>6,"VII"=>7,
					"VIII"=>8,"IX"=>9,"X"=>10,"XI"=>11,"XII"=>12,"XIII"=>13,
					"XIV"=>14,"XV"=>15,"XVI"=>16);
if($chr_name_switch{$chr_test})
{
	$chromosome="chr".$chr_name_switch{$chr_test};
}
print $chromosome,"\n";
my %cassette_elements=("chr15"=>"721500-723000","chr5"=>"116000-117000");


#$chromosome=~s/chr//;
#$chromosome="chr".$chr_name_switch{$chromosome};

sub length_CIGAR
{
	my $CIGARstring=$_[0];
	my $read_line=$_[1];
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
			print $read_line,"\n";
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

sub is_from_cassette
{
	my ($CIGAR,$read_start,$donor_start,$donor_end)=@_;
	my $is_from_cassette=0;
	if($CIGAR=~/^\d+S/)
	{
		if(abs($read_start-$donor_start)<=6)
		{
			$is_from_cassette=1;
		}
	}
	elsif($CIGAR=~/\d+S$/)
	{
		if(abs($read_start+&length_CIGAR($CIGAR)-$donor_end)<=6)
		{
			$is_from_cassette=1;
		}
	}
	return $is_from_cassette;
}

print "START searching for construct reads...\n";
open(INSAM,$input_sam);
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
		my $seqpair_name=$read_column[0];
		my $FLAG=$read_column[1];
		my $seq_self_chr=$read_column[2];
		my $seq_self_map_start=$read_column[3];
		my $seq_CIGAR=$read_column[5];
		my $seq_mate_chr=$read_column[6];
		my $seq_mate_map_start=$read_column[7];
		my $fragment_size=abs($read_column[8]);
		my $this_seq_overlap_target=FALSE;
		my $mate_seq_not_same_location=FALSE;
		#if the read overlapps target interval, returns $this_seq_overlap_target=TRUE
		if($seq_self_chr eq $chromosome)
		{
			if(($donor_start>=$seq_self_map_start)&&($donor_start<($seq_self_map_start+&length_CIGAR($seq_CIGAR,$read_line)))||(($donor_end>=$seq_self_map_start)&&($donor_end<($seq_self_map_start+&length_CIGAR($seq_CIGAR))))||(($donor_start<$seq_self_map_start)&&($donor_end>=($seq_self_map_start+&length_CIGAR($seq_CIGAR)))))
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
