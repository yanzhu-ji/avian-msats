#!/usr/bin/perl
#
# Description: this script takes in a fasta file with flanking regions, and concatenate the flanking regions into one "artificial" sequence. Note that input file is in specific order: flanks of the same target are written one after another.
#
# Usage: $0 <input_flanks.fasta> > <output_concatenated.fasta>
#
# YJ, yji@fieldmusuem.org, 9.27.18


use strict;

my $n = 1;
my ( $id_1, $id_2, $pos_1, $pos_2, $seq_1, $seq_2 );

while ( my $line = <> ){
    chomp $line;
    if ( $n == 1 ){
        # parse seqID
        $line =~ /(.*)_([0-9]{1,}_[0-9]{1,})$/;
        $id_1 = $1;
        $pos_1 = $2;
    }elsif ( $n == 2 ){
        # save seq
        $seq_1 = $line; 
    }elsif ( $n == 3 ){
        # parse seqID 
        $line =~ /(.*)_([0-9]{1,}_[0-9]{1,})$/;
        $id_2 = $1;
        $pos_2 = $2;
    }else {
        #save seq, print new seq and reset $n
        $seq_2 = $line;
        if ( $id_1 eq $id_2 ){
            print STDOUT "$id_1"."_$pos_1"."_$pos_2\n";
            print STDOUT "$seq_1"."$seq_2\n";
            $n = 1;
            next;
        }
    }
    $n++ ; 
}

