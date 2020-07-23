#!/usr/bin/perl

use Bio::AlignIO;
use strict;


($#ARGV == 1) or die <<HERE

    Takes MSA from stdin and produces new MSA in different format to stdout.
  
    Usage: $0 input-format output-format

    Accepted formats: 'clustalw', 'emboss', 'metafasta','fasta', 'maf',
       'mega', 'meme', 'msf', 'nexus', 'pfam', 'phylip', 'psi', 
       'stockholm', 'selex' and 'xmfa'.

HERE
    ;


my $inputformat = shift;
my $outputformat = shift;


my $in  = Bio::AlignIO->new(-fh   => \*STDIN,
			 -format => $inputformat);
my $out = Bio::AlignIO->new(-fh   => \*STDOUT ,
			 -format => $outputformat);

while ( my $aln = $in->next_aln() ) {
    $out->write_aln($aln);
}

