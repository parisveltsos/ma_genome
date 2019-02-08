#!/usr/bin/perl

#Translate a nucleotide seauence to protein
#Input: Fasta file
#Output: Translated sequences


use warnings;
use strict;

use Bio::SeqIO;

my $in_file = $ARGV[0];
die "Enter a fasta file\n" if !$in_file;

my $sequences = Bio::SeqIO->new( 
    -file   => $in_file,
    -format => "fasta",
);


while ( my $dna = $sequences->next_seq ){
    my $protein = $dna->translate( 
        -codontable_id => 1, # standard genetic code
        -frame         => 0, #reading-frame offset 0
    );
    print ">".$dna->display_id, "\n";
    print $protein->seq, "\n";
}
