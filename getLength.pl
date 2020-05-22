#!/usr/bin/perl

# get length information for each sequence in a fasta file. Outputs a two-column file with sequence ID and sequence length.

# Usage: $0 <fasta_file> ( > <output_file> )

use strict;
use Bio::DB::Fasta;

my $fastaFile = shift;
my @ids;
# get @ids by sequence order.  -10.10.19
open (FASTA, "<", $fastaFile) or die ("Cannot open fasta file: $fastaFile!");
while ( my $line = <FASTA> ){
    if ( $line =~ /^>/ ){
        $line =~ s/^>//;
        chomp $line;
        push @ids, $line;
    }
}
close FASTA;

my $db = Bio::DB::Fasta->new( $fastaFile );
#my @ids = $db->ids;
my $n;

foreach my $id ( @ids ){
#	my $length = $db->length( $id );
#	print STDERR "$id\n";
        my $seq = $db->seq( $id );
	$seq =~ s/-//g;
	$seq =~ s/n//g; # added 10/21/2015, to calculate length from consensus
	my $length = length ( $seq );
	my $header = $db-> header ( $id );
	if ( !defined $length ){
		print STDERR "Warning: $id not found! Weird!\n";
		next;
	}
	print STDOUT "$header\t$length\n";  #not sure why header was printed instead of ID...but it did work.
#	print STDERR "$seq\n";
	$n++;
}
print STDERR "Number of seqs in $fastaFile: $n\n";
