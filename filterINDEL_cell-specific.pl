#!usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

sub show_version{
    warn <<END;
    Version 1.1.0 (copyright R. van Boxtel 2017)
END
    exit;
}

sub usage{
    warn <<END;

    Usage:
    Run by typing: perl filterINDEL_cell-specific.pl -i [input VCF] -s [column test sample] -c [column control sample] -out [output directory]

    Required parameters:
    -i | input VCF (including the PATH)
    -s | column number, which contains the genotype field of the TEST sample in the VCF (0-based)
    -c | column number, which contains the genotype field of the CONTROL sample in the VCF (0-based)
    -out | output directory

    Optional parameters:
    -q | minimal GATK quality score for filtering (Default: 250)
    -RD | minimal read depth at a variant position in both TEST and CONTROL sample to pass the filter (Default: 20)
    -VAF | minimal variant allele frequency in the TEST sample for a variant to pass the filter (Default: 0.3)
    -GQ | sample-specific genotype quality as determined by GATK (Default: 99)

END
    exit;
}

die usage() if @ARGV == 0;

my $vcf;
my $sample;
my $control;
my $out_dir;
my $qual = 250; #optional parameter
my $min_rd = 20;
my $min_vaf = 0.1; #optional parameter
my $GQ = 99; #optional parameter
my $help;
my $version;
GetOptions(
    'i=s' => \$vcf,
    's=s' => \$sample,
    'c=s' => \$control,
    'out=s' => \$out_dir,
    'q=s' => \$qual,
    'RD=s' => \$min_rd,
    'VAF=s' => \$min_vaf,
    'GQ=s' => \$GQ,
    'help' => \$help,
    'version' => \$version,
);

if ($help){
    usage();
}

if ($version){
    show_version();
}


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
	next if ($data[5] < $qual); #minimal quality requirement
	next if ($data[6] ne "PASS"); #only consider passed variants

	my @info_control = split(":", $data[$control]); #get specs for control sample
	next if ($info_control[0] ne "0/0"); #remove calls in ref sample
#	next if ($info_control[3] < $GQ); #remove calls with low GQ-score in control sample
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
	next if ($RD_control < $min_rd); #remove positions with less than 20 informative reads in ref sample
	next if ($AA_control > 0); #remove positions with any evidence in ref sample

	my @info_sample = split(":", $data[$sample]); #get specs for test sample
	next if ($info_sample[0] eq "0/0" or $info_sample[0] eq "./."); #remove lines w/o call in test sample
	next if ($info_sample[3] < $GQ); #remove calls with low GQ-score in test sample
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
	next if ($RD_sample < $min_rd); #remove positions with less than 20 informative reads in test sample
	my $AA = $RD_sample - $alleles_sample[0];
	my $VAF = $AA/$RD_sample;

	next if ($AA == 0); #remove calls with no alternative reads, but a call in sample
	next if ($VAF < $min_vaf); #remove samples with low VAF

	print OUT $line, "\n";
    }
}
close IN;
close OUT;