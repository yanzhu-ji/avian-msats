#!/usr/bin/perl
#
# this script sorts the trf .dat file to make sure the starting positions and ending positions of each line is properly sorted.
# Usage: $0 <trf_output.dat> > <trf_output_sorted.dat>
#
# 4/26/2017
#
use strict;

my $seq_line;
my (@entries, @sorted_entries);

while ( my $line = <> ) {
  chomp $line;
  if ( $line =~ /^Sequence/ ){
    if ( $seq_line ne "" ){
      ## sort previously stored lines and print
      @sorted_entries = sort {
        $a->{'start'} <=> $b->{'start'} ||
        $a->{'end'} <=> $b->{'end'}
      } @entries;

      print "$seq_line\n";
      for my $entry ( @sorted_entries ){
        print "$entry->{line}\n";
      }

      ## empty storage
      @entries=();
      @sorted_entries=();
      $seq_line = "";  
    }
      ## save seqID
      $seq_line = $line;
  }else {
    my @field = split " ", $line;
    if ( scalar @field == 15 ) {  
      ## store line
      push @entries, {
        'start' => $field[0],
        'end' => $field[1],
        'line' => $line
      };
    }
  }
}

## sort lastly stored lines and print
print "$seq_line\n";
@sorted_entries = sort {
  $a->{'start'} <=> $b->{'start'} ||
  $a->{'end'} <=> $b->{'end'}
} @entries;

for my $entry ( @sorted_entries ){
  print "$entry->{line}\n";
}
