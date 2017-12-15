#!/usr/bin/perl -w
use strict;
use Vcf;
my $VCF1 = shift;
my $VCF2 = shift;

#open IN1,"gzip -dc $VCF1|" or die $!;
#open IN2,"gzip -dc $VCF2|" or die $!;
#my $line1=<IN1>;
#my $line2=<IN2>;
my $vcfSample = Vcf->new( file => $VCF1, version => '4.1' );
$vcfSample->parse_header();
my (@samples) = $vcfSample->get_samples();
$vcfSample->close();

my $vcfResult = Vcf->new( file => $VCF1, version => '4.1' );
my $vcfRef    = Vcf->new( file => $VCF2, version => '4.1' );

$vcfResult->parse_header();
$vcfRef->parse_header();

my $x = $vcfResult->next_data_hash();
my $y = $vcfRef->next_data_hash();

my %ori;
my %total_heter;
my %switch_error;
my %heter_cnt;
my %wrong_geno;
foreach my $sampleID (@samples) {
	$ori{$sampleID}          = 2;
	$total_heter{$sampleID}  = 0;
	$switch_error{$sampleID} = 0;
	$heter_cnt{$sampleID}    = 0;
	$wrong_geno{$sampleID}   = 0;
}
while ( defined($x) and defined($y) ) {

#	print "process $$x{POS} and $$y{POS}\n";
	if ( $$x{POS} < $$y{POS} ) {
		$x = $vcfResult->next_data_hash();
		next;
	}
	elsif ( $$x{POS} > $$y{POS} ) {
		$y = $vcfRef->next_data_hash();
		next;
	}

	foreach my $sampleID (@samples) {

		#	$vcfResult->set_samples(include=>[$sampleID]);
		#	$vcfRef->set_samples(include=>[$sampleID]);

		my ( $al1, $asep, $al2 ) = return_alleles($vcfResult,$x, $sampleID );

		my ( $bl1, $bsep, $bl2 ) = return_alleles($vcfRef,$y, $sampleID );

#		print "sent:",$al1,$al2,":",$bl1,$bl2,"ori",$ori{$sampleID},"\n";

		compare_geno( $al1, $al2, $bl1, $bl2, $sampleID );
	}
	$x = $vcfResult->next_data_hash();
	$y = $vcfRef->next_data_hash();
}
$vcfResult->close();
$vcfRef->close();

foreach my $sampleID (@samples) {
	print "SampleID:", $sampleID, "\tTotoal heter:", $heter_cnt{$sampleID},
	  "\ttotal switch error:", $switch_error{$sampleID},
	  "\trate:", $switch_error{$sampleID} / ( $heter_cnt{$sampleID} - 1 ),
	  "\ttotal geno error:", $wrong_geno{$sampleID},"\n";
}
sub return_alleles
{
    my ($self,$rec,$column) = @_;
    if ( !exists($$rec{gtypes}) || !exists($$rec{gtypes}{$column}) ) { $self->throw("The column not present: '$column'\n"); }

    my $gtype = $$rec{gtypes}{$column}{GT};
    if ( !($gtype=~$$self{regex_gt}) ) { $self->throw("Could not parse gtype string [$gtype] [$$rec{CHROM}:$$rec{POS}]\n"); }
    my $al1 = $1;
    my $sep = $2;
    my $al2 = $3;
    return ($al1,$sep,$al2);
}
sub compare_geno {
	my ( $alleleA, $alleleB, $alleleC, $alleleD, $sampleID ) = @_;
#	print "received:",$alleleA, $alleleB, $alleleC, $alleleD,"ori",$ori{$sampleID},"\n";
	if ( ( $alleleA + $alleleB ) != ( $alleleC + $alleleD ) ) {
		$wrong_geno{$sampleID}++;
		return 0;
	}
	else {

		if ( ( $alleleC + $alleleD ) == 1 ) {
			$heter_cnt{$sampleID}++;
			if ( $ori{$sampleID} == 2 ) {
				if   ( $alleleC == $alleleA ) { $ori{$sampleID} = 1; }
				else                          { $ori{$sampleID} = 0; }
			}
			if ( ($ori{$sampleID} == 1 &&  $alleleA == $alleleC )
				|| ( $ori{$sampleID} == 0 && $alleleC == $alleleB ) )
			{
				return;
			}
			else {
#				print "error:",$alleleA, $alleleB, $alleleC, $alleleD,"ori",$ori{$sampleID},"\n";
				$switch_error{$sampleID}++;
				$ori{$sampleID} = 1 - $ori{$sampleID};
				return;
			}
		}
	}
}
