#!usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

#-get options
my $file;
my $sample;
my $control;
my $out_dir;
GetOptions(
    'i=s' => \$file,
    's=s' => \$sample,
    'c=s' => \$control,
    'out=s' => \$out_dir,
) or die "ERROR: parameters missing (-i <INPUT>, -s <SAMPLE COLUMN>, -c <CONTROL COLUMN> or -our <OUTPUT DIRECTORY>)\n";

my @types = split("_", $file);
my $type = $types[-1];

my $out = 0;

#-get sample names
open (IN, $file);
while (my $line = <IN>){
    chomp $line;
    next if ($line =~ m/\#\#/);
    if ($line =~ m/\#/){
	my @data = split("\t", $line);
	my @samples = split("_", $data[$sample]);
	my @controls = split("_", $data[$control]);
	$out = $samples[0]."_".$controls[0]."_".$type;
    }
}

close IN;

#-generate output files
open (IN, $file);
open (OUT, '>'.$out_dir."/".$out);

my $stats = $out;
$stats =~ s/.vcf/_stats.txt/;

open (STATS, '>'.$out_dir."/".$stats);
print STATS "CHROM\tSTART\t\tEND\t\tSIZE\tREAD-PAIRS\t\tJUNCTIONS\t\tSPLITS SUM\tPAIRS SUM\n";

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
	my $chr2 = $info[5];
	my $size = $pos2 - $data[1];
	next if ($mapq < 60); #sufficient quality

	my @info_control = split(":", $data[$control]);
	next if ($info_control[3] ne "PASS"); #position must be correctly surveyed in control
	next if ($info_control[0] ne "0/0"); #filter calls in control
	next if ($info_control[9] != 0 or $info_control[11] != 0); #filter any evidence in control

	my @info_sample = split(":", $data[$sample]);
	next if ($info_sample[0] eq "0/0" or $info_sample[0] eq "./."); #only calls in sample
	next if ($info_sample[3] ne "PASS");

	my $pairs = 0;
	my $splits = 0;
	if ($info_sample[8] + $info_sample[9] != 0){
	    $pairs = $info_sample[9] / ($info_sample[8] + $info_sample[9]);
	}
	if ($info_sample[10] + $info_sample[11] != 0){
	    $splits = $info_sample[11] / ($info_sample[10] + $info_sample[11]);
	}

        my $splits_sum = $info_sample[10] + $info_sample[11];
        my $pairs_sum = $info_sample[8] + $info_sample[9];

	if (($splits >= 0.3) or ($pairs >= 0.3 and $size >= 500)) {
	    print STATS $data[0], "\t", $data[1], "\t", $pos2, "\t", $size, "\t", $splits, "\t", $pairs, "\t", $splits_sum, "\t\t", $pairs_sum, "\n";
	    print OUT $line, "\n";
	}
    }
}

close IN;
close OUT;
close STATS;