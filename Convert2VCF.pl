#!/usr/bin/perl -w
my $hap = shift;
my $outVCF = shift;
open IN, "<$hap" or die $!;
print "##fileformat=VCFv4.2
###FILTER=<ID=PASS,Description=\"All filters passed\">
###bcftoolsVersion=1.4.1+htslib-1.4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
my $sampleSize = (split /\s+/,<IN>)[1];
for(my $i=1;$i<=$sampleSize/2;++$i)
{
	print "\t",$i
}
print "\n";
while(<IN>)
{
	last if /segsites/;
}
my $pos_line=<IN>;
chomp($pos_line);
my @sites=split /\s/,$pos_line;
@sites = @sites[1..$#sites];
my @phasedHaps;
my @AFs=(0)x@sites;
while(my $line=<IN>)
{
	chomp($line);
	my @allelesA=split "",$line;
	$line=<IN>;
	chomp($line);
	my @allelesB=split "",$line;
	my @localHaps;
	for(my $i=0;$i<=$#allelesA;++$i)
	{
		push @localHaps,$allelesA[$i]."|".$allelesB[$i];
		$AFs[$i]+=($allelesA[$i]+$allelesB[$i]);
	}
	push @phasedHaps,\@localHaps;
}
close IN;

for(my $i = 0; $i <=$#sites;++$i)
{
	print "chr1\t",int($sites[$i]*@sites*10000),"\t",$i,"\tA\tT\t100\tPASS\tAF=",$AFs[$i]/$sampleSize,"\tGT";
	for(my $j = 0;$j <$sampleSize/2;++$j)
	{
		print "\t",$phasedHaps[$j][$i];
	}
	print "\n";
}

