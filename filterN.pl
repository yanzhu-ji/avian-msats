#!/usr/bin/perl

# Description: filter sequences with number of Ns larger than input in fasta sequence; better use countN.pl to check firstly.
#
# USAGE: $0 -s|sequence	[fasta] 
# 	    -n|number	[integer]
# YJ.10.29.2014, modified from countN.pl

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::PrimarySeq;

my ( $fasta, $n_max, $outfile );

GetOptions (
	's|sequence=s' => \$fasta,
	'n|number=i' => \$n_max,
);

if ( $fasta eq "" || $n_max eq "" ){
	die "Please check the usage!\n";
}

$fasta =~ /(.*)\.fa.*/;

$outfile = $1."_filtered-N.fasta";

my @length_info;
my $inseq = Bio::SeqIO->new(-file => $fasta, 
			    -format => 'Fasta',
			   );
my $outseq = Bio::SeqIO->new(-file => ">$outfile",
			     -format => 'Fasta',
			   );
while (my $seq = $inseq->next_seq() ){
	my $seq_string = $seq->seq();
#d	my $seq_length = $seq->length();
	my $N_count = $seq_string =~ tr/N/N/;
	if ( $N_count <= $n_max ){
		$outseq->write_seq($seq);
	}

#d	my $not_N = $seq_length - $N_count;
#d	push @length_info, [$seq_length, $N_count, $not_N];
#	push @seq_lengths, $seq_length;
#	push @not_Ns, $not_N;
#	print "element done\n";
}
#my $list_N = join (',', @N_counts);
#print "count of Ns:\n$list_N\n\n";

#my $list_seq = join (',', @seq_lengths );
#print "Seq lengths:\n$list_seq\n\n";

#my $list_nN = join (',', @not_Ns );
#print "Remaining:\n$list_nN\n";
#d print "Seq.length\tN.length\tnonN.length\n";
#t foreach my $i (@length_info){
#t	print "$i->[0]\t$i->[1]\t$i->[2]\n";	
#	print "$seq_lengths[$i]\t$N_counts[$i]\t$not_Ns[$i]\n";
#t}

