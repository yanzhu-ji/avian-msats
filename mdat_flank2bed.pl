#!/usr/bin/perl
#
# Description: this script takes an input of the format of mdat (modified dat 
# from Tandem Repeat Finder), and outputs a bed-formatted file with flanking 
# regions of a given length.
#
# Usage: $0     -i|input* <xx.mdat>
#               -l|length* <integer>
#               -t|table* <tab-delimited file with ref. seq lengths>
#
# Yanzhu Ji yanzhuji20@gmail.com 9.11.2017

use strict;
use Getopt::Long;

my ( $infile, $flanking_length, $length_table );

GetOptions(
    "i|input=s" => \$infile,
    "l|length=i" => \$flanking_length,
    "t|table=s" => \$length_table,
);

if ( $infile eq "" || $flanking_length eq "" || $length_table eq "" ){
    die ("Please check the usage!\n");
}

open (TABLE, "<", $length_table) or die ("Cannot open table with sequence length$!\n");

my %hseq_length;
while ( my $line = <TABLE> ){
    chomp $line;
    my @length_fields = split "\t", $line;
    $hseq_length{ $length_fields[0] } = $length_fields[1];
}

#foreach  my $id ( keys %hseq_length ){
#    print "$id -> $hseq_length{$id} \n";
#}

open (INFILE, "<", $infile ) or die ("Cannot open input mdat file:$!\n");

while ( my $line = <INFILE> ){
    chomp $line;
    my @fields = split "\t", $line;
#    print "$fields[0]\n";

    my $start = $fields[1];
    my $end = $fields[2];

    if ( $start >= $flanking_length + 1){
        my $flank_5_start = $start - $flanking_length - 1;
        my $flank_5_end = $start - 1;
        print STDOUT "$fields[0]\t$flank_5_start\t$flank_5_end\n";
    }

    my $seq_length = $hseq_length{ $fields[0] };
    #print "length: $seq_length\n";

    if ( $end <= $seq_length - $flanking_length ){
        my $flank_3_start = $end;
        my $flank_3_end = $end + $flanking_length;
        print STDOUT "$fields[0]\t$flank_3_start\t$flank_3_end\n";
    }
}

