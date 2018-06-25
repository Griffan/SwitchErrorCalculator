#!/usr/bin/perl -w
use strict;
my $VCF1=shift;
my $VCF2=shift;
open IN1,"gzip -dc $VCF1|" or die $!;
open IN2,"gzip -dc $VCF2|" or die $!;
my $line1=<IN1>;
my $line2=<IN2>;
my $ori=2;
my $total_heter=0;
my $switch_error=0;
my $heter_cnt=0;
while(defined($line1) and defined($line2))
{
	if($line1=~/^#/)
	{
		$line1=<IN1>;
		next;
	}
	if($line2=~/^#/)
	{
		$line2=<IN2>;
		next;
	}
	#$line1=<IN1>;
	#$line2=<IN2>;
	chomp $line1;
	chomp $line2;
	my @a=split /\s+/,$line1;
	my @b=split /\s+/,$line2;
	if($a[1] < $b[1])
	{
		$line1=<IN1>;
		next;
	}
	if($a[1] > $b[1])
	{
		$line2=<IN2>;
		next;
	}

	$line1=<IN1>;
	$line2=<IN2>;
	my @A=split /\:/,$a[$#a];
	my @B=split /\:/,$b[$#b];
	$a[$#a]=~/^(.)(\||\/)(.)/;
	my $alleleA=$1;
	my $alleleB=$3;
	if($alleleA ne $alleleB) {$heter_cnt+=1;}
	$b[$#b]=~/^(.)(\||\/)(.)/;
	my $alleleC=$1;
	my $alleleD=$3;
#	print $alleleA,$alleleB,$alleleC,$alleleD,"\n";
	if(($alleleA+$alleleB) != ($alleleC+$alleleD))
	{
		print $a[1],"\t",$a[7],"\t",$a[$#a],"\t",$b[$#b],"\tNO\n";
	}
	else
	{
		 print $a[1],"\t",$a[7],"\t",$a[$#a],"\t",$b[$#b],"\tYES";
		 if(($alleleA+$alleleB) ==1)
		 {
			 $total_heter++;
			 if($ori==2)
			 {
			 	 if($alleleC==$alleleA) {$ori=1;}
			 	 else {$ori=0;}
				 
			 }
			 if(($ori==1 && $alleleA==$alleleC)||($ori==0 && $alleleC==$alleleB))
			 {
				print "\tori:",$ori,"\n";
				 next;
			 }
			 else
			 {
				print "\t Inconsistent $ori $alleleA $alleleB $alleleC\n";
				 $switch_error++;
				 $ori=1-$ori;
				 next;
			 }
		 }
		 else
		 {
		 print "\n";
	 	 }
	}

}
close IN1;
close IN2;
print "Totoal heter:",$heter_cnt,"\ttotal switch error:",$switch_error,"\trate:",$switch_error/($heter_cnt-1),"\n";
