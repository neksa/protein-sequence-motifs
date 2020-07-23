#!/usr/bin/perl

# Author: Simon Mitternacht 2009

# See sub usage for usage instructions.

# --- REMARKS ---

# Found no documentation for MSF-format: parsing assumes that MSA
# starts after a line beginning with '//' and that each block of
# sequence data starts with at least one empty line, followed by a
# line with two numbers (contents not checked, line is simply
# ignored). Each line of the sequence data consists of one column with
# sequence-identifer and then 1-5 columns of sequences. Gaps are
# assumed to be indicated by '.'. Columns are assumed to be separated
# by white-space ('\s+').
#

use Getopt::Long;
use strict;

sub usage { 
    my $basename =`basename $0`;
    chomp $basename;
    print <<"HERE";

    Calculates a frequency matrix from a multiple sequence alignment
    (MSA) in MSF-format.  Either the profile of the whole MSA can be
    used, or the environment of a conserved residue. In the latter
    case all gaps are removed from the analysis.

    Input is taken from stdin and output written to stdout.

    Usage: $basename [--help] [--whole_msa] [--conserved c --before b --after a]

HERE
    exit;
}

### handle command line arguments ###
($#ARGV > -1) or usage; # there always has to be at least on argument

my $display_help = 0;

my $interval_defined = 0;
my $before = -1;
my $after = -1;
my $start = 0;
my $end = 0;
my $conserved = -1;
my $use_whole_msa = 0;
my $length = 0;

GetOptions("conserved=i"=>\$conserved,
	   "before=i"=>\$before,
	   "after=i"=>\$after,
	   "whole_msa!"=>\$use_whole_msa,
	   "help!"=>\$display_help
	   );
if ($display_help) {
    usage;
}

if ($conserved != -1 || $before != -1 || $after != -1) {
    $interval_defined = 1;
}
if ($conserved != -1 && $before != -1 && $after != -1) {
    $interval_defined = 2;
    --$conserved; # numbering in MSF starts from 1, internally it starts from 0
}
if ($interval_defined && $use_whole_msa) {
    die "Error: cannot use both whole MSA and an interval, make a choice\n";
}
if ($interval_defined == 1) {
    die "Error: sequence interval not correctly defined.\n";
}


# Get the MSA (could use Bioperl::AlignIO here instead..., not
# installed on all BCCS computers?)
my @msa = get_sequences_from_msf(*STDIN); 

# gaps are included to allow consistency check
my $residues = "ACDEFGHIKLMNPQRSTVWY."; 

#### print and calculate output ####
my $pos = 0;
my $N_seq = $#msa+1;

# use whole alignment if nothing is specified
if ($interval_defined == 0) {
    $start = 0;
    $end = length $msa[0]; # assume all sequences have the same length
    $length = $end + 1;
}
# otherwise use the specified segment
else {
    $start = $conserved - $before;
    $end = $conserved + $after;
    $length = $before + $after + 1;
    my @tmp_msa;
    foreach (@msa) {
	push @tmp_msa, get_seq_surrounding_align_pos($_, $conserved, $start, $end);
    }
    @msa = @tmp_msa;
}

# print header
print "PROFILE 3\n";
print "BEGIN\n";
print "MATRIX K=$N_seq L=$length\n";

printf("%4d",$length);
for (my $i = 0; $i < 21; ++$i) {
    print "    ".(substr($residues, $i, 1));
}
print "\n";

# output profile
for (my $i = 0; $i < $length; ++$i) {
    my %profile;
    # add up frequencies
    foreach (@msa) {
	my $res = substr($_, $i, 1);
	++$profile{$res};
    }
    printf ("%4d", $pos);
    # print frequencies
    for (my $j = 0; $j < 21; ++$j) {
	my $res = substr($residues, $j, 1);
	printf("%5d", $profile{$res});
    }
    print "\n";
    ++$pos;
}


print "END\n";

# get the sequence surrounding a certain conserved position (i.e. remove the gaps)
sub get_seq_surrounding_align_pos
{
    # $seq should be a string, the other three arguments numbers
    my ($seq, $conserved, $start, $end) = @_;
    my @seq_array = split '', $seq;
    my @result;
    my $offset = $conserved-$start;
    my $length = $end-$start+1;
    ($conserved > $start) or
	die "Error: conserved residue has to be within interval.\n\!";

    $result[$offset] = $seq_array[$conserved];
    my $pos = $conserved-1;
    
    # check that the conserved position does not have a gap
    unless ( $result[$offset] =~ m/[A-Za-z]/) { 
	die "Error: conserved residue has gap.\n!";
    }

    # fill up before ...
    for (my $i = $offset-1; $i >= 0 && $pos >= 0;) {
	my $res = $seq_array[$pos];
	# if not gap insert position into seq
	if ($res =~ m/[A-Za-z]/) {
	    $result[$i--] = $res;
	} 
	--$pos;
    }
    # ... and after conserved position
    $pos = $conserved+1;
    for (my $i = $offset+1; $i < $length && $pos <= $#seq_array;) {
	my $res = $seq_array[$pos];
	if ($res =~ m/[A-Za-z]/) {
	    $result[$i++] = $res;
	} 
	++$pos;
    }
    return join '', @result;
}

#### extract aligned sequences from file ####
sub get_sequences_from_msf
{
    my $input = shift;
    # use hash while parsing
    my %sequences;
    # return array in the end
    my @result;

    # flags
    my $started = 0;
    my $newblock = 0;

    while(<$input>) {
	# find start of MSA
	if ($started == 0) {
	    if (m/^\/\// == 1) {
		$started = 1;
	    }
	    next;
	}    
	
	# every sequence block starts with an empty line ...
	if (m/^\s*$/) {
	    $newblock = 1; 
	    next;
	}
	# ... and a line with sequence numbering
	if ($newblock == 1) {
	    $newblock = 0;
	    next;
	}
	
	# remove leading white-space
	chomp;
	$_ =~ s/^\s+//;
	my @columns = split /\s+/;
	
	#append sequence
	for (my $i = 1; $i <= $#columns; ++ $i) {
	    $sequences{$columns[0]} .= $columns[$i];
	}
    }
    # copy sequences to array
    foreach (keys %sequences) {
	push @result, $sequences{$_};
    }
    return @result;
}


