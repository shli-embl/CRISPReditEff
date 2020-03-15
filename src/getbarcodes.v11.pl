#usage: perl getbarcode.v2.pl input.sam output.tbl startpos seqlength righttrim lefttrim 
use strict;
#@env is a global variable
my @env;
sub Node
{
	my $_node = { A => undef, T => undef, C => undef, G => undef, key => undef, string => undef };
#set multiple value at sam0e time
	my %params = @_;
	map { $_node->{$_} = $params{$_} if defined $params{$_} } keys %$_node;
	
	return $_node;

}

sub getAllPath
{
	my $tmp=$_[0];

	if(!$tmp)
	{ 
		return;
	}
	if(($tmp->{"key"})&&(!$tmp->{"A"})&&(!$tmp->{"T"})&&(!$tmp->{"C"})&&(!$tmp->{"G"}))
	{
		push(@env, $tmp->{"string"});
		
	}
#This optional command output all PATHs in the tree and their frequency
#	print OUT $tmp->{"string"},"\t",$tmp->{"key"},"\n";
	getAllPath($tmp->{"A"});
	getAllPath($tmp->{"T"});
	getAllPath($tmp->{"C"});
	getAllPath($tmp->{"G"});
}

sub getTScore
{
	my $tmp=$_[0];
	my $string=$_[1];
	my $score=$tmp->{"A"}->{"key"}+$tmp->{"T"}->{"key"}+$tmp->{"C"}->{"key"}+$tmp->{"G"}->{"key"};
	for (my $i=0;$i<length($string);$i++)
	{
		my $offspring=0;
		my @nuc=qw{A T C G};
		foreach my $nucl(@nuc)
		{
			if($tmp->{$nucl})
			{
				$offspring++;

			}

		}
		if($offspring>1)
		{
			foreach my $nucl2(@nuc)
			{
				if(($nucl2 ne substr($string,$i,1))&&($tmp->{$nucl2}))
				{
					$score=$score-$tmp->{$nucl2}->{"key"};
				}
			}
			
		}
		$tmp=$tmp->{substr($string,$i,1)};
	}
	return $score;
}

sub getUScore
{
	my $tmp=$_[0];
	my $string=$_[1];
	my $score=0;
	for (my $i=0;$i<length($string);$i++)
	{
		my $offspring=0;
		my @nuc=qw{A T C G};
		foreach my $nucl(@nuc)
		{
			if($tmp->{$nucl})
			{
				$offspring++;

			}

		}
		if(($offspring>1)||(!$score))
		{
			
			$score=$tmp->{substr($string,$i,1)}->{"key"};
		}
		$tmp=$tmp->{substr($string,$i,1)};
	}
	return $score;
}

sub aln_mm
{
	my $seq1=$_[0];
	my $seq2=$_[1];
	# split the sequences
	my $match=10;
	my $mismatch=-10;
	my $gop=-10;
	my $gep=-1;
	my @smat;
	my @tb;

	#evaluate substitutions
	my $len1=length($seq1);
	my $len2=length($seq2);
	#first base must match
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
			my $start_penalty=0;
			if(($i>5)&&($j>5))
			{
				$start_penalty=0;
			}
			#calcul du score
			if(substr($seq1,$i-1,1) eq substr($seq2,$j-1,1))
			{
				$s=$match;
			}
			elsif((substr($seq1,$i-1,1) eq "N")||(substr($seq2,$j-1,1) eq "N"))
			{
				$s=0;
			}
			else
			{
				$s=$mismatch+$start_penalty;
			}
		
			my $sub=$smat[$i-1][$j-1]+$s;
			my $del;
			if($i==$len1)
			{
				$del=$smat[$i  ][$j-1]+$start_penalty;
			}
			else
			{
				if($tb[$i][$j-1]==-1)
				{
					$del=$smat[$i  ][$j-1]+$gep+$start_penalty;
				}
				else
				{
					$del=$smat[$i  ][$j-1]+$gop+$gep+$start_penalty;
				}
			}
			
			my $ins;
			if($tb[$i-1][$j]==1)
			{
				$ins=$smat[$i-1][$j]+$gep+$start_penalty;
			}
			else
			{
				$ins=$smat[$i-1][$j]+$gop+$gep+$start_penalty;
			}
			
			
			
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
	#Output en Fasta:
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
			if($aln1[$i+1] eq "-")
			{}
			else
			{
				#gap in first 5 bases discared
				if(($aln_len-1-$i)>5)
				{
					$mm++;
				}
				else
				{
					$mm=1000;
				}
			}
		}
		elsif($aln2[$i] eq "-")
		{
			if($aln2[$i+1] eq "-")
			{}
			else
			{
				if(($aln_len-1-$i)>5)
				{
					$mm++;
				}
				else
				{
					$mm=1000;
				}
			}
		}

	}
#	print "\n",$seq1,"\n";
#	print "mismatch:",$mm,"\ttraceback:",$smat[$len1][$len2],"\n";
#	for(my $k=$aln_len-1;$k>=0;$k--)
#	{
#		print $aln1[$k];
#	}
#	print "\n";
#	for(my $k=$aln_len-1;$k>=0;$k--)
#	{
#		print $aln2[$k];
#	}
#	print "\n";
	return $mm;
}


#define the tree root;
my $rootleft = Node();
my $rootright = Node();

#process the SAM file to get all substrings
open(IN,$ARGV[0]);
open(OUT,">$ARGV[1]");
#eg. 3994
my $start=$ARGV[2];
#eg. 31
my $length=$ARGV[3];
while(my $content=<IN>)
{
	if($content!~/^\@/)
	{
		my @tmp=split("\t",$content);
		my $flag=$tmp[1];
		my $pos=$tmp[3];
		my $alncode=$tmp[5];
		my $seq=$tmp[9];
#filter flag
		if(($flag eq "83")||($flag eq "99")||($flag eq "147")||($flag eq "163"))
		{
#overlap the left
			my $left=0;
			if($alncode=~/^(\d+)S.+/)
			{
				$left=$1;
			}
			my $readlength=0;
			my $tmpread=$alncode;
			while($tmpread)
			{
				if($tmpread=~/^(\d+)(\w)/)
				{
					my $c=$1;
					my $sym=$2;
					if($sym ne "D")
					{
						$readlength+=$c;
					}
					$tmpread=~s/^\d+\w//;
				}
				else
				{
					print "There is error ALNCODE\n";
					exit 0;
				}
				
			}
			if((($pos-$left)>=($start-$readlength))&&(($pos-$left)<($start+$length)))
			{

				my @read;
				my $pread=0;
				my $pgenome=$pos-1;
				while($alncode=~/^(\d+)(\w)/)
				{
					if($2 eq "S")
					{
						for(my $i=0;$i<$1;$i++)
						{
							$read[$pread]=$pgenome+0.5;
							$pread++;
						}
					}
					elsif($2 eq "M")
					{
						for(my $i=0;$i<$1;$i++)
						{
							
							$read[$pread]=$pgenome+1;
							$pread++;
							$pgenome++;
						}
					}
					elsif($2 eq "D")
					{
						for(my $i=0;$i<$1;$i++)
						{
							$pgenome++;
						}
					}
					elsif($2 eq "I")
					{
						for(my $i=0;$i<$1;$i++)
						{
							
							$read[$pread]=$pgenome+0.5;
							$pread++;
							
						}
					}
					
					$alncode=~s/^\d+\w//;
				}
#left and right go tree align
#left	

				if($pos <= $start)
				{
					my $tmp=$rootleft;
					my $i=0;
					while(($read[$i]<=($start-1))&&($i<=$readlength))
					{
						
						$i++;
					}
#process the sequence					
					my $j=0;
					my $string=substr($seq,$i,$readlength-$i);
#					print $content,"\n",$i,"\t",$readlength,"\t",$read[$i],"\t",$read[$i+1],"\t",$string,"\n\n";
					while($j<length($string))
					{

						if($tmp->{substr($string,$j,1)})
						{
							my $tmpstring=$tmp->{"string"};

							$tmp=$tmp->{substr($string,$j,1)};
							$tmp->{"key"}++;
							$tmp->{"string"}=$tmpstring.substr($string,$j,1);
						}
						else
						{
							my $tmpstring=$tmp->{"string"};
							$tmp->{substr($string,$j,1)}=Node();
							$tmp=$tmp->{substr($string,$j,1)};
							$tmp->{"key"}=1;
							$tmp->{"string"}=$tmpstring.substr($string,$j,1);
						}
						$j++;

					}
					
				}
				else
				{
#The reads overlaped with right side as well
					my $tmp=$rootright;
					my $i=$readlength-1;
					while(($read[$i]>=($start+$length))&&($i>=0))
					{
						
						$i--;
					}
#process the sequence					

					my $j=$i;
					my $string=substr($seq,0,$i+1);
					while($j>=0)
					{

						if($tmp->{substr($string,$j,1)})
						{
							my $tmpstring=$tmp->{"string"};

							$tmp=$tmp->{substr($string,$j,1)};
							$tmp->{"key"}++;
							$tmp->{"string"}=$tmpstring.substr($string,$j,1);
						}
						else
						{
							my $tmpstring=$tmp->{"string"};
							$tmp->{substr($string,$j,1)}=Node();
							$tmp=$tmp->{substr($string,$j,1)};
							$tmp->{"key"}=1;
							$tmp->{"string"}=$tmpstring.substr($string,$j,1);
						}
						$j--;
					}
				}
			}
		}
	}
}





undef @env;
#get all paths in the left tree, stored in @env
#get scores for paths
my %lefttrim;
my %righttrim;
my %suspeciousTleft;
my %suspeciousUleft;
getAllPath($rootleft);
my %pathsTleft;
my %pathsUleft;


#my $righttrim="GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTTTTTTGTTTTTTATGTCTTTTTGAAGAGTGC";
my $righttrim=$ARGV[4];
#my $lefttrim="CACTACGAAATCTTGAGATCGGGCGTTCGACTCGCCCCCGGGAGAGATGGCCGGCATGGTCCCAGCCTCCTCGCTGGCGCCGGCTGGGCAACACCTTCGGGTGGCGAATGGGACTTTGGCGCGCC";
my $lefttrim=$ARGV[5];
foreach my $x (@env)
{


	if(length($x) > $length)
	{	
		my $flag=1;
		my $pt=13;
		while(($pt<length($x))&&$flag)
		{
#			print "What:","\n";
			if(&aln_mm(substr($x,$pt-1,length($x)-$pt+1),$righttrim) <= 0.1*(length($x)-$pt+1))
			{
				$flag=0;
				$pathsTleft{substr($x,0,$pt-1)}=getTScore($rootleft,substr($x,0,$pt-1));
				$pathsUleft{substr($x,0,$pt-1)}=getUScore($rootleft,substr($x,0,$pt-1));
				$righttrim{substr($x,0,$pt-1)}="right_trim";
#				print "test:",substr($x,0,$pt-1),"\t",$x,"\n";
			}
			else
			{
				$pt++;
			}
		
		}
		if($flag)
		{
			$suspeciousTleft{substr($x,0,$pt-1)}=getTScore($rootleft,substr($x,0,$pt-1));
			$suspeciousUleft{substr($x,0,$pt-1)}=getUScore($rootleft,substr($x,0,$pt-1));
		}
#		print substr($x,0,$length),"\t",getTScore($rootleft,substr($x,0,$length)),"\t",getUScore($rootleft,substr($x,0,$length)),"\n";
	}
	else
	{
		my $newx=$x;
		while(length($newx)<$length)
		{
			$newx=$newx."-";
		}
		$pathsTleft{$newx}=getTScore($rootleft,$x);
		$pathsUleft{$newx}=getUScore($rootleft,$x);

	}
}


#get all paths in the right tree, stored in @env
#get scores for paths
undef @env;
getAllPath($rootright);
my %pathsTright;
my %pathsUright;
my %suspeciousTright;
my %suspeciousUright;
foreach my $x(@env)
{


	if(length($x) > $length)
	{	

		my $flag=1;
		my $pt=13;
		while(($pt<length($x))&&$flag)
		{
			
			if(&aln_mm(substr($x,$pt-1,length($x)-$pt+1),(my $rvlt=reverse($lefttrim))) <= 0.1*(length($x)-$pt+1))
			{
				$flag=0;
				$pathsTright{my $rvst=reverse(substr($x,0,$pt-1))}=getTScore($rootright,substr($x,0,$pt-1));
				$pathsUright{my $rvst=reverse(substr($x,0,$pt-1))}=getUScore($rootright,substr($x,0,$pt-1));
				$lefttrim{my $rvst=reverse(substr($x,0,$pt-1))}="left_trim";
	#			print "left:",substr($x,0,$pt-1),"\t",$x,"\n";
			}
			else
			{
				$pt++;
			}
		
		}
		if($flag)
		{
			$suspeciousTright{my $rvst=reverse(substr($x,0,$pt-1))}=getTScore($rootright,substr($x,0,$pt-1));
			$suspeciousUright{my $rvst=reverse(substr($x,0,$pt-1))}=getUScore($rootright,substr($x,0,$pt-1));
		}
			
		
	#	$pathsTright{reverse(substr($x,0,$length))}=getTScore($rootright,substr($x,0,$length));
	#	$pathsUright{reverse(substr($x,0,$length))}=getUScore($rootright,substr($x,0,$length));
	}
	else
	{
		my $newx=reverse($x);
		while(length($newx)<$length)
		{
			$newx="-".$newx;
		}
		$pathsTright{$newx}=getTScore($rootright,$x);
		$pathsUright{$newx}=getUScore($rootright,$x);

	}
}


#merge left and right tree
my %mergeTpaths;
my %mergeUpaths;
foreach (keys %pathsTleft)
{
	$mergeTpaths{$_}+=$pathsTleft{$_};
	$mergeUpaths{$_}+=$pathsUleft{$_};
}
foreach (keys %pathsTright)
{
	$mergeTpaths{$_}+=$pathsTright{$_};
	$mergeUpaths{$_}+=$pathsUright{$_};
}

#foreach my $x(keys %mergeTpaths)
#{
	#foreach my $y(keys %mergeTpaths)
	#{
	#	if($x ne $y)
	#	{
	#		my $flag=1;
#			my $n=0;
	#		my $match=0;
		#	while($n<$length)
	#		{
		#		print $n,"\t",substr($x,$n,1),"\t",substr($y,$n,1),"\n";
		#		if((substr($x,$n,1) ne substr($y,$n,1))&&(substr($x,$n,1) ne "-")&&(substr($y,$n,1) ne "-"))
		#		{
		#			$flag=0;
		#		}
	#			if((substr($x,$n,1) eq substr($y,$n,1))&&(substr($x,$n,1) ne "-")&&(substr($y,$n,1) ne "-"))
	#			{
	#				$match++;
	#			}
	#			$n++;
				
	#		}
	#		if($flag&&($match>=3))
	#		{

	#			my $combinestr;
	#			for(my $i=0;$i<$length;$i++)
	#			{
	#				if(substr($x,$i,1) eq "-")
	#				{
		##				$combinestr=$combinestr.substr($y,$i,1);
		#			}
	#				else
	#				{
	#					$combinestr=$combinestr.substr($x,$i,1);
	#				}
	#			}
				
	#			$mergeTpaths{$combinestr}=$mergeTpaths{$x}+$mergeTpaths{$y};
	#			$mergeUpaths{$combinestr}=$mergeUpaths{$x}+$mergeUpaths{$y};
	#			delete $mergeTpaths{$y};
	#			delete $mergeTpaths{$x};
	#			delete $mergeUpaths{$y};
	#			delete $mergeUpaths{$x};
	#			
	#		}
#		}
#	}
#}
#print OUT %mergeTpaths,"\n";
#print OUT %mergeUpaths,"\n";
my $num=1;
#foreach my $key (sort {$mergeTpaths{$b}<=>$mergeTpaths{$a}} keys %mergeTpaths)
foreach my $key (sort {($mergeTpaths{$b}<=>$mergeTpaths{$a}) || ($mergeUpaths{$b}<=>$mergeUpaths{$a})} keys %mergeTpaths)
{
	if(length($key)<($length+20))
	{
		print OUT "SEQ$num:\t",$key,"\tTscore:",$mergeTpaths{$key},"\tUscore:",$mergeUpaths{$key},"\t",$righttrim{$key},"\t",$lefttrim{$key},"\n";
		$num++;
	}
	
	
}
print OUT "\n\nSuspecious sequences:\n";
foreach my $key (sort {($mergeTpaths{$b}<=>$mergeTpaths{$a}) || ($mergeUpaths{$b}<=>$mergeUpaths{$a})} keys %mergeTpaths)
{
	if(length($key)>=($length+20))
	{
		print OUT "SEQ$num:\t",$key,"\tTscore:",$mergeTpaths{$key},"\tUscore:",$mergeUpaths{$key},"\t",$righttrim{$key},"\t",$lefttrim{$key},"\n";
		$num++;
	}
	
#	
}
print OUT "\nleft:\n";
foreach my $key (sort {($suspeciousTleft{$b}<=>$suspeciousTleft{$a}) || ($suspeciousUleft{$b}<=>$suspeciousUleft{$a})} keys %suspeciousTleft)
{

	print OUT "SEQ$num:\t",$key,"\tTscore:",$suspeciousTleft{$key},"\tUscore:",$suspeciousUleft{$key},"\n";
	
	$num++;
}
print OUT "\nright:\n";
foreach my $key (sort {($suspeciousTright{$b}<=>$suspeciousTright{$a}) || ($suspeciousUright{$b}<=>$suspeciousUright{$a})} keys %suspeciousTright)
{

	print OUT "SEQ$num:\t",$key,"\tTscore:",$suspeciousTright{$key},"\tUscore:",$suspeciousUright{$key},"\n";
	
	$num++;
}

#	

undef @env;
undef $rootleft;
undef $rootright;

