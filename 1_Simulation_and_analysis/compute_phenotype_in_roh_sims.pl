#!/usr/bin/perl -w

use strict;
use POSIX qw/floor/;
use POSIX qw/ceil/;

if($#ARGV < 3){
    print STDERR "<muts file> <roh file> <vcf file> <out file>" and die;
}


my $FEATUREFILE = $ARGV[0];
my $ROHFILE = $ARGV[1];
my $VCFFILE = $ARGV[2];
my $OUTFILE = $ARGV[3];
my @vcfFileList;

push(@vcfFileList,$VCFFILE);

my %effect;
my %h;
my $REF = "REF";
$effect{$REF} = 0;
$h{$REF} = 0;
open(FIN,"gunzip -c $FEATUREFILE |") or die $!;
while(my $line = <FIN>){
	chomp $line;
	my ($pos, $type, $mid, $h, $z, $s, $junk) = split(/\s+/,$line,7);
	$effect{$mid} = $z;
	$h{$mid} = $h;	
}
close(FIN);

my %ROH; #popind -> chr -> array of arrays
  
print STDERR "Reading $ROHFILE...\n";
open(ROH,"gunzip -c $ROHFILE |") or die $!;
   	
my $pop;
my $ind;
#my %ind2pop;
while (my $rohline = <ROH>){
	chomp $rohline;

	if($rohline =~ m/^track .+Ind: (.+) Pop:(.+) ROH.+/){
		$ind = $1;
		$pop = $2;
#		$ind2pop{$ind} = $pop;
	}
	else{
		my ($chr, $start, $end, $class, $size, $junk) = split(/\s+/,$rohline,6);
		my @tmp = ($start,$end-1,$class);
		push(@{$ROH{$ind}{$chr}},\@tmp);
	}
}
close(ROH);
  	

my %stabilizing;
my %positive;
#load individuals, all the same so just take chr1
open(FIN,"gunzip -c $vcfFileList[0] |") or die $!;
print STDERR "Loading individual list.\n";
my @indlist;
my ($junk1, $junk2, $junk3, $junk4, $junk5, $junk6, $junk7, $junk8, $junk9);
for my $line (<FIN>){
	chomp $line;
	if($line =~ m/^##/){
	    next;
	}
	elsif($line =~ m/^#CHROM/){
	    ($junk1, $junk2, $junk3, $junk4, $junk5, $junk6, $junk7, $junk8, $junk9, @indlist) = split(/\s+/,$line);

	    for my $ind (@indlist){
			$stabilizing{$ind}{'A'} = 0;
			$stabilizing{$ind}{'B'} = 0;
			$stabilizing{$ind}{'C'} = 0;
			$stabilizing{$ind}{'NONE'} = 0;
			$positive{$ind}{'A'} = 0;
			$positive{$ind}{'B'} = 0;
			$positive{$ind}{'C'} = 0;
			$positive{$ind}{'NONE'} = 0;
	    }
	}
	else{
	    last;
	}
}
close(FIN);

for (my $chr = 1; $chr <= 1; $chr++){
	my $chrstr = "chr$chr";
	my $file = $vcfFileList[$chr-1];
	open(FIN,"gunzip -c $file |") or die $!;
	print STDERR "${chrstr}\n";
	my $count = 0;
	while(my $line = <FIN>){
		chomp $line;
		if($line =~ m/^#/){
		    next;
		}

		my ($c, $pos, $rsid, $mid, $q, $filter, $info, $format, @genotypes) = split(/\s+/,$line);
		
		if($count % 100000 == 0){
			print STDERR "$chrstr $pos\n";
		}
		$count++;

		my @mids = split(/,/,$mid);
		my %allele;
		$allele{0} = $REF;
		my $a = 1;
		for my $m (@mids){
			$allele{$a} = $m;
			$a+=1;
		}

		if(exists $effect{$mid}){
			#print STDERR "$chr $pos ", join(' ',@genotypes), "\n";
		    for (my $i = 0; $i <= $#indlist; $i++){
				my ($a1, $a2) = split(/\|/,$genotypes[$i]);
				if($a1 eq '.'){
			    	next;
				}
				if($a1 eq $a2)
				{
			    	my $class = hitsInterval(\@{$ROH{$indlist[$i]}{$chrstr}},$pos);
			    	my $effect = $effect{$allele{$a1}};
			    
				    if($class eq '0'){
						$stabilizing{$indlist[$i]}{'NONE'} += $effect;
						$positive{$indlist[$i]}{'NONE'} += abs($effect);
				    }
			    	else{
				    	$stabilizing{$indlist[$i]}{$class} += $effect;
				    	$positive{$indlist[$i]}{$class} += abs($effect);
			    	}
				}
				if($a1 ne $a2)
				{
			    	my $class = hitsInterval(\@{$ROH{$indlist[$i]}{$chrstr}},$pos);
				    #print STDERR "$class\n";
				    my $effect = $effect{$allele{$a1}};
				    my $h = $h{$allele{$a1}};
				    if($class eq '0'){
						$stabilizing{$indlist[$i]}{'NONE'} += $h*$effect;
						$positive{$indlist[$i]}{'NONE'} += $h*abs($effect);
				    }
				    else{
					    $stabilizing{$indlist[$i]}{$class} += $h*$effect;
					    $positive{$indlist[$i]}{$class} += $h*abs($effect);
				    }
				    $effect = $effect{$allele{$a2}};
				    $h = $h{$allele{$a2}};
					    if($class eq '0'){
						$stabilizing{$indlist[$i]}{'NONE'} += $h*$effect;
						$positive{$indlist[$i]}{'NONE'} += $h*abs($effect);
				    }
				    else{
					    $stabilizing{$indlist[$i]}{$class} += $h*$effect;
					    $positive{$indlist[$i]}{$class} += $h*abs($effect);
				    }
				}
			}
		}
	}
	close(FIN);
}

open(FOUT,">",$OUTFILE) or die $!;

print FOUT "phenoAstab phenoBstab phenoCstab phenoNONEstab phenoApos phenoBpos phenoCpos phenoNONEpos\n";

for my $ind (@indlist){
		print FOUT "$ind ";
		print FOUT $stabilizing{$ind}{'A'}, " ";
		print FOUT $stabilizing{$ind}{'B'}, " ";
		print FOUT $stabilizing{$ind}{'C'}, " ";
		print FOUT $stabilizing{$ind}{'NONE'}, " ";
		print FOUT $positive{$ind}{'A'}, " ";
		print FOUT $positive{$ind}{'B'}, " ";
		print FOUT $positive{$ind}{'C'}, " ";
		print FOUT $positive{$ind}{'NONE'}, "\n";
}


close(FOUT);


#takes an array of intervals and a query
#returns the 3rd element of the interval if query is inside one of the intervals
#returns 0 otherwise
sub hitsInterval 
  {
    my @intervals = @{$_[0]};
    my $query = $_[1];
    #print $padding, "\n";
    my $numIntervals = @intervals;

    if($numIntervals == 0)
      {
	return 0;
      }

    #print STDERR "number of intervals: $numIntervals\n";
    
    #print "Starting search at ", floor($numIntervals/2), "\n";

    my @range = (0,$numIntervals-1);
    my $index = floor(($range[1]-$range[0])/2);
    my $result;
    my $inside = 0;

    do
      {
	#print "$query\n";
	#print "Index: $index\n";
	#print "Range: ", $range[0], " ", $range[1], "\n";
       
	if ($query <= $intervals[$index][1] && #$query is inside an interval
	    $query >= $intervals[$index][0])
	  {
	    #print STDERR "$query is inside [",$intervals[$index][0], ",",$intervals[$index][1] ,"]\n";
	    $result = 0;
	    $inside = 1;
	  }
	elsif ($query > $intervals[$index][1])
	  {
	    #print STDERR "$query > ", $intervals[$index][1], "\n";
	    if ($index == $numIntervals-1)
	      {
		$result = 0;
	      }
	    elsif($query < $intervals[$index+1][0])
	      {
		$result = 0;
		#print STDERR "Off target query: ";
		#print STDERR "[",$intervals[$index][0], ",",$intervals[$index][1] ,"]";
		#print STDERR " < $query < ";
		#print STDERR "[",$intervals[$index+1][0], ",",$intervals[$index+1][1] ,"]\n";
	      }
	    else
	      {
		$result = 1;
		$range[0] = $index;
		$index = ceil(($range[1]-$range[0])/2)+$range[0];
		#print "New index = ", $index, "\n";
		#print "[",$intervals[$index][0], ",",$intervals[$index][1] ,"]\n";# and die;
	      }
	  }
	elsif ($query < $intervals[$index][0])
	  {
	    #print STDERR "$query < ", $intervals[$index][0], "\n";
	    if ($index == 0)
	      {
		$result = 0;
	      }
	    elsif($query > $intervals[$index-1][1])
	      {
		$result = 0;
		#print STDERR "Off target query: ";
		#print STDERR "[",$intervals[$index-1][0], ",",$intervals[$index-1][1] ,"]";
		#print STDERR " < $query < ";
		#print STDERR "[",$intervals[$index][0], ",",$intervals[$index][1] ,"]\n";
	      }
	    else
	      {
		$result = -1;
		$range[1] = $index;
		$index = floor(($range[1]-$range[0])/2)+$range[0];
	      }
	  }
    } while ($result != 0);

    #print "$result\n";
    #die;
    #print STDERR $result, "\n";

    if($inside == 0)
      {
	return $inside;
      }
    
    return $intervals[$index][2];
    
  }
