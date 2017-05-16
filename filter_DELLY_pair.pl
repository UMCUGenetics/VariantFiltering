#!usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

sub show_version{
    warn <<END;
    Version 1.0.0 (copyright R. van Boxtel 2017)
END
    exit;
}

sub usage{
    warn <<END;

    Usage:
    Run by typing: perl filter_DELLY_pair.pl -i [input VCF] -s [column test sample] -c [column control sample] -out [output directory]

    Required parameters:
    -i | input VCF (including the PATH)
    -s | column number, which contains the genotype field of the TEST sample in the VCF (0-based)
    -c | column number, which contains the genotype field of the CONTROL sample in the VCF (0-based)
    -out | output directory

    Optional parameters:
    -MQ | minimal MAPQ score (Default: 60)
    -GQ | minimal read depth at a variant position in both TEST and CONTROL sample to pass the filter (Default: 99)

END
    exit;
}

die usage() if @ARGV == 0;

#-get options
my $file;
my $sample;
my $control;
my $out_dir;
my $min_mapq = 60; #optional parameter
my $min_gq = 99; #optional parameter
my $version;
my $help;
GetOptions(
    'i=s' => \$file,
    's=s' => \$sample,
    'c=s' => \$control,
    'out=s' => \$out_dir,
    'MQ=s' => \$min_mapq,
    'GQ=s' => \$min_gq,
    'version' => \$version,
    'help' => \$help,
) or die "ERROR: parameters missing (-i <INPUT>, -s <SAMPLE COLUMN>, -c <CONTROL COLUMN> or -our <OUTPUT DIRECTORY>)\n";

if ($help){
    usage();
}

if ($version){
    show_version();
}

#--get names op the sample and type of SV
my @types = split("_", $file);
my $type = $types[-1];
$type =~ s/.vcf//;

my $sample_name = "Error";
my $control_name = "Error";

open (IN, $file);
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
my $out = $sample_name."_".$control_name."_".$type."_noEvidenceControl.vcf";

#-filtering
open (IN, $file);
open (OUT, '>'.$out_dir."/".$out);

my $stats = $out;
$stats =~ s/.vcf/_stats.txt/;

open (STATS, '>'.$out_dir."/".$stats);
print STATS "CHROM\tSTART\t\tEND\t\tSIZE\tTYPE\tCHROM2\tVARIANT_PAIRS\tTOTAL_PAIRS\tVARIANT_SPLITS\tTOTAL_SPLITS\n";

#-filter vcf
while (my $line = <IN>){
    chomp $line;
    if ($line =~ m/\#/){
	print OUT $line, "\n";
    }
    else{
	my @data = split("\t", $line);
	next if ($data[6] ne "PASS"); #line must be passed
	my @info = split(";", $data[7]);
	my $pos2 = 0;
	my $mapq = 0;

	foreach my $element (@info){
	    if ($element =~ m/END/){
		my @ends = split("=", $element);
		$pos2 = $ends[1];
	    }
	    elsif ($element =~ m/MAPQ/){
		my @mapqs = split("=", $element);
		$mapq = $mapqs[1];
	    }
	}
	my (undef, $chr2) = split('=\s*', $info[5]);
	my $size = $pos2 - $data[1];
	next if ($mapq < $min_mapq); #sufficient quality

	my @info_control = split(":", $data[$control]);
#	next if ($info_control[3] ne "PASS"); #position must be correctly surveyed in control
	next if ($info_control[0] ne "0/0"); #filter calls in control
	next if ($info_control[2] < $min_gq); #filter if genotype quality is smaller than defined
	next if ($info_control[9] != 0 or $info_control[11] != 0); #filter any evidence in control

	my @info_sample = split(":", $data[$sample]);
	next if ($info_sample[0] eq "0/0" or $info_sample[0] eq "./."); #only calls in sample
	next if ($info_sample[2] < $min_gq); #filter if genotype quality is smaller than defined
#	next if ($info_sample[3] ne "PASS");

        my $pairs_sum = $info_sample[8] + $info_sample[9];
        my $splits_sum = $info_sample[10] + $info_sample[11];

	print STATS $data[0], "\t", $data[1], "\t", $pos2, "\t", $size, "\t", $type, "\t", $chr2, "\t", $info_sample[9], "\t", $pairs_sum, "\t", $info_sample[11], "\t", $splits_sum, "\n";
	print OUT $line, "\n";

    }
}

close IN;
close OUT;
close STATS;