use strict;
use Getopt::Long;
#@env is a global variable
my @env;
my $display_help;
my $input_bam;
my $output_file;
#coordinate of the sequence for detemination on reference (e.g. coord of guide start)
my $position;
#right adjacent sequence (for identification of boundary); recommended length equals to maximum read length. 
my $right_trim;
#left adjacent sequence
my $left_trim;

#load parameters
GetOptions('h' => \$display_help, 'in=s' => \$input_bam, 'out=s' =>\$output_file, 'pos=s' =>\$position, 'R=s' =>\$right_trim, 'L=s' =>\$left_trim);
if($display_help)
{
	print "Command: \n\tperl extractGuideBarcode.pl [-h] -in INPUT_BAM -out OUTPUT_BAM -chr CONTIG_NAME -pos CHR:START_COORD-END_COORD -R RIGHT_FLANKING_SEQ -L LEFT_FLANKING_SEQ\n\n";
	print "Function: \n\tIdentifying Guide/Barcode sequence from mapped reads (bams). The sequence on reference is represented as N strings.\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput BAM.\n";
	print "\t-out\tOutput result table.\n";
	print "\t-pos\tContig name and position of the variable sequence. (e.g. REDI_contig:1483-1502)\n";
	print "\t-R\tFlanking sequence to the right of variable sequence, will be used for determination of boundary. An recommended length of this sequence should equal to the maximum read length in the provided bam.\n";
	print "\t-L\tFlanking sequence to the left of variable sequence, will be used for determination of boundary. An recommended length of this sequence should equal to the maximum read length in the provided bam.\n";
	exit 1;
}

if((!$input_bam)||(!$output_file)||(!$position)||(!$right_trim)||(!$left_trim))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

if($position!~/.+\:\d+\-\d+$/)
{
	print "Not correct format for -pos argument. (contig_name:start_pos-end_pos)\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}
my @split_pos=split(":",$position);
my $contig_name=$split_pos[0];
my @split_coord=split("-",$split_pos[1]);

#define subfunctions
#Define new node
sub Node
{
	my $_node = { A => undef, T => undef, C => undef, G => undef, key => undef, string => undef };
#set multiple value at same time
	my %params = @_;
	map { $_node->{$_} = $params{$_} if defined $params{$_} } keys %$_node;
	
	return $_node;

}

#Traverse the tree; record path (sequence) in given array
sub getAllPath
{
	#$string_path is an array reference
	my ($tmp_node,$ref_to_list)=@_;

	if(!$tmp_node)
	{ 
		return;
	}
	if(($tmp_node->{"key"})&&(!$tmp_node->{"A"})&&(!$tmp_node->{"T"})&&(!$tmp_node->{"C"})&&(!$tmp_node->{"G"}))
	{
		$ref_to_list->{"list"}->{$tmp_node->{"string"}}++;
		$ref_to_list->{"test"}="try that";
	}
	
#print OUT $$tmp_node->{"string"},"\t",$$tmp_node->{"key"},"\n";
#This optional command output all PATHs in the tree and their frequency
	getAllPath($tmp_node->{"A"},$ref_to_list);
	getAllPath($tmp_node->{"T"},$ref_to_list);
	getAllPath($tmp_node->{"C"},$ref_to_list);
	getAllPath($tmp_node->{"G"},$ref_to_list);
}

sub getTScore
#Given a startnode and the represented sequence of it, calculate the T (Total) score of it, or total number of visits
{
	my ($tmp_node,$string)=@_;
	#initialize score=sum of reads passing all sub-nodes of root (sum of A,T,C,G)
	my $score=$tmp_node->{"A"}->{"key"}+$tmp_node->{"T"}->{"key"}+$tmp_node->{"C"}->{"key"}+$tmp_node->{"G"}->{"key"};
	#iteration: moving a node downward every time
	for (my $i=0;$i<length($string);$i++)
	{
		my $offspring=0;
		my @nuc=qw{A T C G};
		foreach my $nucl(@nuc)
		{
			if($tmp_node->{$nucl})
			{
				$offspring++;
			}
		}
		#if a parent node has only one offspring node, score remains unchanged (parent visits = offspring visits); else if a parent node has >1 offspring nodes, score subtracts the number of visits from offsprings visits with unmatched bases (other three nucleotides except the one included in the given path)
		if($offspring>1)
		{
			foreach my $nucl2(@nuc)
			{
				if(($nucl2 ne substr($string,$i,1))&&($tmp_node->{$nucl2}))
				{
					$score=$score-$tmp_node->{$nucl2}->{"key"};
				}
			}			
		}
		$tmp_node=$tmp_node->{substr($string,$i,1)};
	}
	return $score;
}

sub getUScore
#Given a startnode and the represented sequence of it, calculate the U (Unique) score of it, or unique number of visits
#The outcome is the No. of visits to an query path on its last diversification from other paths
{
	my ($tmp_node,$string)=@_;
	#initialize the score as 0;
	my $score=0;
	#iteration: moving a node downward every time
	for (my $i=0;$i<length($string);$i++)
	{
		my $offspring=0;
		my @nuc=qw{A T C G};
		foreach my $nucl(@nuc)
		{
			if($tmp_node->{$nucl})
			{
				$offspring++;
			}
		}
		#if is the root or a parent node with >1 offspring nodes, update score to the No. of visits at the query path
		if(($offspring>1)||(!$score))
		{			
			$score=$tmp_node->{substr($string,$i,1)}->{"key"};
		}
		$tmp_node=$tmp_node->{substr($string,$i,1)};
	}
	return $score;
}

#semi-global NW alignemtn; aln_mm(a,b): a is sequence of test, b is the flanking sequence; partially match for b is not penalized
sub aln_mm
{
	my ($seq1,$seq2)=@_;
	my $match=10;
	my $mismatch=-10;
	my $gop=-10;
	my $gep=-1;
	#score matrix
	my @smat;
	#trace back matrix
	my @tb;
	my $len1=length($seq1);
	my $len2=length($seq2);
	#initialization of score matrix
	$tb[0][0]=-10;
	for (my $i=0; $i<=$len1; $i++)
	{
		$smat[$i][0]=0;
		$tb[$i][0 ]= 1;
	}
	for (my $j=0; $j<=$len2; $j++)
	{
		$smat[0][$j]=0;
		$tb[0 ][$j]=-1;
	}
	for (my $i=1; $i<=$len1; $i++)
	{
		for (my $j=1; $j<=$len2; $j++)
		{
			my $s;
			if(substr($seq1,$i-1,1) eq substr($seq2,$j-1,1))
			{
				$s=$match;
			}
			elsif((substr($seq1,$i-1,1) eq "N")||(substr($seq2,$j-1,1) eq "N"))
			{
				$s=0;
			}
			
			#Subsititution
			my $sub=$smat[$i-1][$j-1]+$s;
			
			#Deletion
			#CUSTOMIZED:(semi-local alignment) do not penalize gap at end of sequence A
			my $del;
			if($i==$len1)
			{
				$del=$smat[$i][$j-1];
			}
			else
			{
				if($tb[$i][$j-1]==-1)
				{
					$del=$smat[$i][$j-1]+$gep;
				}
				else
				{
					$del=$smat[$i][$j-1]+$gop+$gep;
				}
			}
					
			#Insertion	
			my $ins;
			if($tb[$i-1][$j]==1)
			{
				$ins=$smat[$i-1][$j]+$gep;
			}
			else
			{
				$ins=$smat[$i-1][$j]+$gop+$gep;
			}	
			
			#set trace back: max(sub,del,ins)
			if(($sub>$del) && ($sub>$ins))
			{
				$smat[$i][$j]=$sub;
				$tb[$i][$j]=0;
			}
			elsif($del>$ins)
			{
				$smat[$i][$j]=$del;
				$tb[$i][$j]=-1;
			}
			else 
			{
				$smat[$i][$j]=$ins;
				$tb[$i][$j]=1;
			}
		}
	}
	#Trace back
	my $i=$len1;
	my $j=$len2;
	my $aln_len=0;
	my $mm=0;
	my @aln1;
	my @aln2;
	while (!($i==0 && $j==0))
	{
		if ($tb[$i][$j]==0)
		{
			$aln1[$aln_len]=substr($seq1,$i-1,1);
			$aln2[$aln_len]=substr($seq2,$j-1,1);
			$i--;
			$j--;
		}
		elsif ($tb[$i][$j]==-1)
		{
			$aln1[$aln_len]='-';
			$aln2[$aln_len]=substr($seq2,$j-1,1);
			$j--;			
		}
		elsif ($tb[$i][$j]==1)
		{
			$aln1[$aln_len]=substr($seq1,$i-1,1);
			$aln2[$aln_len]='-';
			$i--;
		}
		$aln_len++;		
	}
	#calculate No. of mismatches
	for ($i=$aln_len-1; $i>=0; $i--)
	{
		if(($aln1[$i] ne "-")&&($aln2[$i] ne "-")&&($aln1[$i] ne $aln2[$i]))
		{
			if(($aln1[$i] ne "N")&&($aln2[$i] ne "N"))
			{
				$mm++;
			}
		}
		elsif($aln1[$i] eq "-")
		{
			if($aln1[$i+1] ne "-")
			{
				#gap in first 5 bases discared; gap regardless of length regarded as one mismatch
				if(($aln_len-1-$i)>5)
				{
					$mm++;
				}
				else
				{
					$mm=10000;
				}
			}
		}
		elsif($aln2[$i] eq "-")
		{
			if($aln2[$i+1] ne "-")
			{
				if(($aln_len-1-$i)>5)
				{
					$mm++;
				}
				else
				{
					$mm=10000;
				}
			}
		}
	}
	return $mm;
}

sub trim_read
{
	
	#trim direction: L: trim on the left, R: trim on the right
	my ($read_seq,$junction_seq,$trim_direction)=@_;
	#left trim
	my $n_base = substr($read_seq,0,length($read_seq));
	my $if_trim = 0;
	my $trimmed_seq = $read_seq;
	if($trim_direction eq "L")
	{
		$read_seq = reverse($read_seq);
		$n_base = reverse($n_base);
		$junction_seq = reverse($junction_seq);
		$trimmed_seq = reverse($trimmed_seq);
	}
	while((!$if_trim)&&($n_base))
	{
		
		#allow 1 error per 10 bases (no gap is allowed in first 5 bases)
		if(&aln_mm($n_base,$junction_seq) <= 0.1 * length($n_base))
		{
			$if_trim = 1;
			$trimmed_seq = substr($read_seq,0,length($read_seq)-length($n_base));
		}
		else
		{
			if(length($n_base) > 1)
			{
				$n_base=substr($n_base,1,length($n_base)-1);
			}
			else
			{
				$n_base="";
			}
		}
	}
	if($trim_direction eq "L")
	{
		$read_seq = reverse($read_seq);
		$n_base = reverse($n_base);
		$junction_seq = reverse($junction_seq);
		$trimmed_seq = reverse($trimmed_seq);
	}
	my $output;
	$output->{"trimmed_seq"}=$trimmed_seq;
	$output->{"if_trimmed"}=$if_trim;
	return $output;
}

sub decimal_to_binary
{
	my $decimal=$_[0];
	my $binary=sprintf("%b",$decimal)	
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

sub add_new_path_to_tree
{
	my ($treeroot,$seq)=@_;
	my $j=0;
	my $tmp_node=$treeroot;
	while($j<length($seq))
	{
		#if an offspring node already exists, move to it and visit++
		if($tmp_node->{substr($seq,$j,1)})
		{
			$tmp_node=$tmp_node->{substr($seq,$j,1)};
			$tmp_node->{"key"}++;
		}
		else
		{
			my $tmp_seq=$tmp_node->{"string"};
			$tmp_node->{substr($seq,$j,1)}=Node();
			$tmp_node=$tmp_node->{substr($seq,$j,1)};
			$tmp_node->{"key"}=1;
			$tmp_node->{"string"}=$tmp_seq.substr($seq,$j,1);
		}
		$j++;
	}
}

sub max
{
	my $tmp;
	foreach my $x (@_)
	{
		if((!$tmp)||($tmp<$x))
		{
			$tmp=$x;
		}
	}
	return $tmp
}

### Here starts the main program ###
#define the tree roots; left-tree starts from the first base of given interval, right-tree starts (reversely) from the last base of the given interval
my $rootleft = Node();
my $rootright = Node();
open(IN,"samtools view -h ".$input_bam." | ");
open(OUT,">".$output_file);
while(my $read_line=<IN>)
{
	if($read_line!~/^\@/)
	{
		my @read_column=split("\t",$read_line);
		my $flag=&decimal_to_binary($read_column[1]);
		my $read_contig=$read_column[2];
		my $read_pos=$read_column[3];
		my $read_CIGAR=$read_column[5];
		my $read_seq=$read_column[9];
		#if Flag has less than 8 binary digits and starts with "11"(not a PCR duplicates, properly mapped reads)
		if((length($flag)<=8)&&(substr(my $rev=reverse($flag),0,2) eq "11")&&($read_contig eq $contig_name))
		{
			my $left_clipped=0;
			my $right_clipped=0;
			my $read_mapped_length=&length_CIGAR($read_CIGAR);
			if($read_CIGAR=~/^(\d+)S.+/)
			{
				$left_clipped=$1;
			}
			if($read_CIGAR=~/(\d+)S$/)
			{
				$right_clipped=$1;
			}
			#if the read start upstream to the tested region & the read spans left junction
			if(($read_pos < ($split_coord[0]+$split_coord[1])/2)&&(($read_pos + $read_mapped_length + $right_clipped) > $split_coord[0]))
			{
				my $ref_left_trim=&trim_read($read_seq,$left_trim,"L");
				#if the left boundary match
				if($ref_left_trim->{"if_trimmed"})
				{
					&add_new_path_to_tree($rootleft,$ref_left_trim->{"trimmed_seq"});
				}
			}
			#if the read start downstream to the tested region & the read spans right junction
			elsif(($read_pos > ($split_coord[0]+$split_coord[1])/2)&&(($read_pos - $left_clipped) < $split_coord[1]))
			{
				my $ref_right_trim=&trim_read($read_seq,$right_trim,"R");
				#if the right boundary match
				if($ref_right_trim->{"if_trimmed"})
				{
					&add_new_path_to_tree($rootright,my $rev =reverse($ref_right_trim->{"trimmed_seq"}));
				}
			}
		}
	}
}
my $combined_call_ref;
my $left_seqs_ref;
$left_seqs_ref->{"list"};
&getAllPath($rootleft,$left_seqs_ref);
foreach(keys %{$left_seqs_ref->{"list"}})
{
	my $ref_right_trimmed_seq=&trim_read($_,$right_trim,"R");
	$combined_call_ref->{$ref_right_trimmed_seq->{"trimmed_seq"}}->{"TScore"}+=&getTScore($rootleft,$ref_right_trimmed_seq->{"trimmed_seq"});
	$combined_call_ref->{$ref_right_trimmed_seq->{"trimmed_seq"}}->{"UScore"}+=&getUScore($rootleft,$ref_right_trimmed_seq->{"trimmed_seq"});
	$combined_call_ref->{$ref_right_trimmed_seq->{"trimmed_seq"}}->{"left_junction"}=1;
	$combined_call_ref->{$ref_right_trimmed_seq->{"trimmed_seq"}}->{"right_junction"}=&max($ref_right_trimmed_seq->{"if_trimmed"},$combined_call_ref->{$ref_right_trimmed_seq->{"trimmed_seq"}}->{"right_junction"});
}

my $right_seqs_ref;
$right_seqs_ref->{"list"};
&getAllPath($rootright,$right_seqs_ref);
foreach(keys %{$right_seqs_ref->{"list"}})
{
	my $ref_left_trimmed_seq=&trim_read(my $rev=reverse($_),$left_trim,"L");
	$combined_call_ref->{$ref_left_trimmed_seq->{"trimmed_seq"}}->{"TScore"}+=&getTScore($rootright,my $rev=reverse($ref_left_trimmed_seq->{"trimmed_seq"}));
	$combined_call_ref->{$ref_left_trimmed_seq->{"trimmed_seq"}}->{"UScore"}+=&getUScore($rootright,my $rev=reverse($ref_left_trimmed_seq->{"trimmed_seq"}));
	$combined_call_ref->{$ref_left_trimmed_seq->{"trimmed_seq"}}->{"left_junction"}=&max($ref_left_trimmed_seq->{"if_trimmed"},$combined_call_ref->{$ref_left_trimmed_seq->{"trimmed_seq"}}->{"left_junction"});
	$combined_call_ref->{$ref_left_trimmed_seq->{"trimmed_seq"}}->{"right_junction"}=1;
}

print OUT "Sequence\t(T)otal_N_of_reads\t(U)nique_N_of_reads\toverlap_5_junction\toverlap_3_junction\n";

foreach(sort {$combined_call_ref->{$b}->{"TScore"}<=>$combined_call_ref->{$a}->{"TScore"}} keys %{$combined_call_ref})
{
	print OUT  $_,"\t",$combined_call_ref->{$_}->{"TScore"},"\t",$combined_call_ref->{$_}->{"UScore"};
	if($combined_call_ref->{$_}->{"left_junction"})
	{
		print OUT "\tY";
	}
	else
	{
		print OUT "\tN";
	}
	if($combined_call_ref->{$_}->{"right_junction"})
	{
		print OUT "\tY";
	}
	else
	{
		print OUT "\tN";
	}
	print OUT "\n";
}
