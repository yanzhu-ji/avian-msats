#!/usr/bin/perl
# Description: this script takes (1) a mdat file (modified dat file that TRF outputs)
#   (2) a tab-delimited file with seq ID and their lengths (3) flanking length and
#   (4) prefix of output files.  The script writes three files that have both, one 
#   and none flanking sequence, given the flanking length.
#   Usage: $0 -i <in.mdat>
#             -t <length.txt>
#             -l <flanking length, integer>
#             -o <output prefix>
#   YJ, yanzhuji20@gmail.com
#   9.13.2017

use strict;
use Getopt::Long;

my ( $infile, $length_table, $flanking_length, $out_prefix );

GetOptions (
    "i|input=s" => \$infile,
    "t|table=s" => \$length_table,
#    "f|flanking_option=s" => \$flanking_option,
    "l|length=i" => \$flanking_length,
    "o|out_prefix=s" => \$out_prefix,
);

open (TABLE, "<", $length_table) or die ("Cannot open table with sequence length$!\n");

my %hseq_length;
while ( my $line = <TABLE> ){
    chomp $line;
    my @length_fields = split "\t", $line;
    $hseq_length{ $length_fields[0] } = $length_fields[1];
}

close TABLE;

my $number = scalar (keys %hseq_length);
print STDERR "Total number of sequences: $number\n";

open (BOTH, ">", $out_prefix."_w-both-flanks.mdat") or die ("Cannot write to output:$!\n");
open (ONE, ">", $out_prefix."_w-one-flank.mdat" ) or die ("Cannot write to output:$!\n");
open (NONE, ">", $out_prefix."_w-none-flank.mdat" ) or die ("Cannot write to output:$!\n");

open (INFILE, "<", $infile) or die  ("Cannot open input mdat file:$!\n");

while ( my $line = <INFILE> ){
    chomp $line;
    my @fields = split "\t", $line;
    my $start = $fields[1];
    my $end = $fields[2];

    my $seq_length = $hseq_length{ $fields[0] };
    if ( $seq_length == "" ) {
        print STDERR "Warning: sequence length for $fields[0] not found!\n";
        next;
    }
    
    if ( $start >= $flanking_length + 1 && $end <= $seq_length - $flanking_length ){
        print BOTH "$line\n";
    } elsif ( $start < $flanking_length + 1 && $end > $seq_length - $flanking_length ){
       print NONE "$line\n";
    } else {
       print ONE "$line\n";
    }
} 

close BOTH;
close ONE;
close NONE;
close INFIlE;

