#!usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

my $vcf;
my $sample;
my $control;
my $out_dir;
GetOptions(
    'i=s' => \$vcf,
    's=s' => \$sample,
    'c=s' => \$control,
    'out=s' => \$out_dir,
) or die "ERROR: parameters missing (-i <INPUT>, -s <SAMPLE_COLUMN>, -c <CONTROL_COLUMN> or -out <OUPUT DIRECTORY>)\n";

#--get names of the samples
my $sample_name = "Error";
my $control_name = "Error";

open (IN, $vcf);
while (my $line = <IN>){
    chomp $line;
    next if ($line =~ m/\#\#/);
    if ($line =~ m/\#/){
	my @data = split("\t", $line);
	$sample_name = $data[$sample];
	$control_name = $data[$control];
    }
}
close IN;

#--generate ASC-specific output file
my $out = $sample_name."_".$control_name."_INDELs_noEvidenceControl.vcf";

#--filtering
open (IN, $vcf);
open (OUT, '>'.$out_dir."/".$out);

while (my $line = <IN>){
    chomp $line;
    my @data = split("\t", $line);
    if ($line =~ m/\#/){
	print OUT $line, "\n";
    }
    else{
	next if ($data[0] eq "X" or $data[0] eq "Y" or $data[0] eq "MT"); #get autosomal genome;
	next if ($data[2] =~ m/rs/ and $data[2] !~ m/COSM/); #remove SNP_ids w/o COSMIC id
	next if ($data[4] =~ m/\,/);
	next if ($data[5] < 100); #minimal quality requirement
	next if ($data[6] ne "PASS"); #only consider passed variants

	my @info_control = split(":", $data[$control]); #get specs for control sample
	next if ($info_control[0] ne "0/0"); #remove calls in ref sample
	my @alleles_control = split(",", $info_control[1]);
	my $RD_control = 0;
	foreach my $allele (@alleles_control){
	    if ($RD_control == 0){
	    $RD_control = $allele;
	    }
	    else{
	    $RD_control = $RD_control + $allele;
	    }
	}
	my $AA_control = $RD_control - $alleles_control[0];
	next if ($RD_control < 20); #remove positions with less than 20 informative reads in ref sample
	next if ($AA_control > 0); #remove positions with any evidence in ref sample

	my @info_sample = split(":", $data[$sample]); #get specs for test sample
	next if ($info_sample[0] eq "0/0" or $info_sample[0] eq "./."); #remove lines w/o call in test sample
	my @alleles_sample = split(",", $info_sample[1]);
	my $RD_sample = 0;
	foreach my $allele (@alleles_sample){
	    if ($RD_sample == 0){
	    $RD_sample = $allele;
	    }
	    else{
	    $RD_sample = $RD_sample + $allele;
	    }
	}
	next if ($RD_sample < 20); #remove positions with less than 20 informative reads in test sample
	my $AA = $RD_sample - $alleles_sample[0];
	my $VAF = $AA/$RD_sample;

	next if ($AA == 0); #remove calls with no alternative reads, but a call in sample
	next if ($VAF < 0.3); #remove samples with low VAF

	print OUT $line, "\n";
    }
}
close IN;
close OUT;