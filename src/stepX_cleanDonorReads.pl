#usage: perl stepX_cleanDonorReads.pl input.sam output.sam donor_start_position donor_length
#donor start: 3654 end:3762 length:109
#Output is a SAM file containing only usable reads (spanning either junction of REDI donor locus) for calling donor sequences from input alignment
use strict;
use constant FALSE => 1==0;
use constant TRUE => not FALSE;
my $input_sam=$ARGV[0];
my $output_sam=$ARGV[1];
my $left_junction=$ARGV[2];
my $right_junction=$left_junction+$ARGV[3];

#a sub function determining the mapped length of a read from CIGAR string
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

#a function determines TRUE or FALSE, whether a read span the junctions
sub span_junction
{
	my ($CIGAR,$read_start,$left_junction,$right_junction)=@_;
	my $span_junction=0;
	my $left_or_right=0;
	my $read_end=$read_start+&length_CIGAR($CIGAR);
	if(($read_start<($left_junction-3))&&($read_end>($left_junction+3)))
	{
		$span_junction=TRUE;
		$left_or_right="L";
	}
	elsif(($read_start<($right_junction-3))&&($read_end>($right_junction+3)))
	{
		$span_junction=TRUE;
		$left_or_right="R";
	}
	if(($read_start<($left_junction-3))&&($read_end>($right_junction+3)))
	{
		$span_junction=TRUE;
		$left_or_right="RT";
	}
	
	return ($span_junction,$left_or_right);
}

#print "START searching for construct reads...\n";
open(INSAM,$input_sam);
my %on_target_list_unique;
my %span_read_count=("L"=>0,"R"=>0,"RT"=>0);
my $line_count=0;
while(my $read_line=<INSAM>)
{
	chomp $read_line;
	if($line_count%100000==1)
	{
#		print "Processed ",$line_count-1," lines.\n";		
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
		my $this_read_span_junction=FALSE;
		my $this_read_left_or_right=0;
		#FLAG < 1024 ~ not a PCR duplicate
		if(($seq_self_chr eq "REDI_cassette")&&($FLAG<1024))
		{
			($this_read_span_junction,$this_read_left_or_right)=&span_junction($seq_CIGAR,$seq_self_map_start,$left_junction,$right_junction);
			
		}
		
		if($this_read_span_junction)
		{
			$on_target_list_unique{$seqpair_name}=1;
			$span_read_count{$this_read_left_or_right}++;

		}
	
	}
}
close INSAM;
#End with a hash table of usable reads


#start filtering input SAM and generate output SAM with only usable reads
#print "Start filtering SAM file\n";
open(INSAM,$input_sam);
open(OUTSAM,">".$output_sam);
$line_count=0;
while(my $read_line=<INSAM>)
{
	chomp $read_line;
	if($line_count%100000==1)
	{
#		print "Processed ",$line_count-1," lines.\n";
	}
	$line_count++;
	if($read_line!~/^\@/)
	{
		my @read_column=split("\t",$read_line);
		my $readpair_name=$read_column[0];
		if($on_target_list_unique{$readpair_name})
		{

			print OUTSAM $read_line,"\n";

		}			
	}
	else
	{
		print OUTSAM $read_line,"\n";
	}
}
print $input_sam,"\t",$span_read_count{"L"}."L".$span_read_count{"R"}."R".$span_read_count{"RT"}."RT","\n";

