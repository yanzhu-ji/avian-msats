#!/usr/bin/perl
# Description: this script can: 
#   (1) option -c/--check: checks for redundant report of tandem repeats reported by Tandem Repeat Finder.
#   (2) option -f/--filter: if redundancy is detected, filter the redundant entries based on alignment score
# Usage: $0 -i <*.dat> 
#           -c/--check, or,
#           -f/--filter 
#
# Yanzhu Ji 12.6.2016

use strict;
use Getopt::Long;

my ( $infile, $check, $filter );
my $seqID;
my ( $this_start, $this_end, $o_start, $o_end, $this_almScore, $o_almScore );
my ( $this_line, $o_line, $best_line);
my ( $best_start, $best_end, $best_almScore);
my ( $outfile, $outfile_uniq );

GetOptions(
    "i|infile=s" => \$infile,
    "c|check" => \$check,
    "f|filter" => \$filter,
);

unless ( $check or $filter ) {
	die ("Please specify at least one of -c and -f!\n" );
}

open (INFILE, "<", $infile ) or die ("Cannot open input file: $!\n");

if ( $filter ){
    $outfile = $infile;
    $outfile_uniq = $infile;
    $outfile =~ s/\.dat/\.mmdat/;
    $outfile_uniq =~ s/\.dat/.mdat/;
    open (OUTFILE, ">", $outfile ) or die ("Cannot create new file: $!\n");
}

while ( $this_line = <INFILE> ) {
    chomp $this_line;

    if ( $this_line =~ /^Sequence/ ) {
        #  compare o_line and this_line and print?
        if ( $best_line ne "" ){
            print OUTFILE "$seqID $best_line\n";
	    #t print "printing line: $best_line\n";
	}

        ## save sequence ID
        my @seqID_field = split " ", $this_line, 3;
        $seqID = $seqID_field[1];
        $best_start = "";
        $best_end = "";
        $best_almScore = "";
        $best_line = ""; 
    } else {
        my @field = split " ", $this_line;
        if ( scalar @field == 15 ) {
            ## save first and second fields, and compare with last memorized fields
            $this_start = $field[0];
            $this_end = $field[1];
            $this_almScore = $field[7];

	    #t print "reading line: $seqID\t$this_start\t$this_end\n";
            #t print "old line: $seqID\t$best_start\t$best_end\n";
            if ( $best_end eq "" ){
                $best_start = $this_start;
                $best_end = $this_end;
                $best_almScore = $this_almScore;
                $best_line = $this_line; 
            } elsif ( $this_start <= $best_end ) { ## if two regions are overlapping
                if ( $check ) {
                    ## report if redundancy is detected
                    print "$seqID\n";
                    print "$best_line\n";
                    print "$this_line\n\n";
                }
                if ( $filter ) {
                    ## compare the alignment scores of two entries and print after reformatting
                    if ( $this_almScore > $best_almScore ) {
                        #t print "Updated scores: $this_almScore\t$best_almScore\n";
                        $best_line = $this_line;
                        $best_start = $this_start;
                        $best_end = $this_end;
                        $best_almScore = $this_almScore;
                        next;
		    }else {
                        #t print "Scores:$this_almScore\t$best_almScore\n";
			next;
			#$best_line = $o_line;
                        #$best_start = $o_start;
                        #$best_end = $o_end;
                        #$best_almScore = $o_almScore;
                    } 
        #            print OUTFILE "$seqID $best_line\n";
                }
            }else { ## if two regions are not overlapping, then print the last result from comparisons
		print OUTFILE "$seqID $best_line\n";
	        #t print "printing line: $best_line\n";
                $best_start = $this_start;
                $best_end = $this_end;
                $best_almScore = $this_almScore;
                $best_line = $this_line;
            }
            #$best_start = $this_start;
            #$best_end = $this_end;
            #$best_almScore = $this_almScore;
            #$best_line = $this_line;
        }
    }
}

print OUTFILE "$seqID $best_line\n";
#t print "printing line: $best_line\n";
close OUTFILE;

if ( $filter ){
    `uniq  $outfile > $outfile_uniq`;
    `rm $outfile`;
}


