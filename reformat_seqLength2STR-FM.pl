#!/usr/bin/perl
## Description: this script is used to parse getLength.pl output into a STR-FM input file. 
##              The getLength.pl output includes two columns: one with fasta ID that includes locus ID and species ID, and a second column of sequence length.
##              This script also takes an integer to limit the locus to print if the number of sequences (lengths) is less than some value. Default: 5
## Usage:       $0  <getLength.pl_output.txt>  <motif-locusID chart>  <interger>   >  <output>

use strict;
use File::Basename;

my $length_info = $ARGV[0];
my $motif_info = $ARGV[1];
my $limit = $ARGV[2];

#  get super_motif
my %super_motif ;
open (SUPERMOTIF, "<", $motif_info ) or die ("Cannot open file:$motif_info!\n");
while ( my $line = <SUPERMOTIF> ){
    chomp $line;
    my @motif_info = split "\t", $line;
    $super_motif{$motif_info[1]} = $motif_info[2];
}
close SUPERMOTIF;

my ($new_qID, $old_qID);
my ($new_motif, $old_motif);
my @lengths = ();

open (INPUT, "<", $length_info ) or die ("Cannot open file: $length_info!\n");
while (my $line = <INPUT> ){
    chomp $line;
    my @fields = split "\t", $line;
    $fields[0] =~ s/-R-//;
    #print STDERR "reading: $fields[0]\n";
    $fields[0] =~ /^([a-zA-Z-]{6,40}-L[0-9]{1,3})-*/;
    $new_qID = $1;
    #print STDERR "reading: $new_qID\n";
    my @temp = split "-", $new_qID;
    my $locus_ID = $temp[-1];
    my $msat_length = $fields[1];
    #push @lengths, $msat_length;
    $new_motif = $super_motif{$locus_ID};

    if ( $old_qID eq "" && $new_qID ne "" ){
        #print "First line found!\n";
        push @lengths, $msat_length;
    }elsif ( $new_qID eq $old_qID ){
        push @lengths, $msat_length;
    }else {
        if ( scalar (@lengths) >= $limit ){
            print "$old_qID\t", join(",",@lengths), "\t$old_motif\n";
            @lengths = ();
            push @lengths, $msat_length;
        }
    }
    $old_qID = $new_qID;
    $old_motif = $new_motif;
}

if ( scalar (@lengths) >= $limit ){
    print "$old_qID\t", join(",",@lengths), "\t$old_motif\n";
}

