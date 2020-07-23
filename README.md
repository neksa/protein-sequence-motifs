# Protein sequence motif derivation with extreme sensitivity

Alexander Goncearenco (developed in 2007-2012 at the University of Bergen, Norway)

Designed to derive sequence motifs of elementary functional loops (EFL) in proteins.

## The package contains:

    * converge - find conserved motifs in sequence
    * clustering - reduce the number of derived motifs by combining the similar ones
    * search - find matches of a predefined motif in sequece database

## Utilities:
    * findPermutatons.R - find positional permutations that destroy relative positional distances and hence the natural motif
    * convertMSA.pl - convert multiple sequence alignment to initial motif
    * composition.csv - average prokaryotic proteomic composition
    * profile_logo.py - generate profile logo graphics
    * profile_analysis.py
    * calc_logo.py - a shorter numpy version of calc logo
    * smooth.py - a PyMOL script to smooth protein chain with a window

Check other branches as well for modifications.

## Changes in this branch
    * removed MPI for searchPSSM and convergePSSM
    * added position slicing from-to
    * added output to variable lenght motif files (format #4)
    * multiple randomizations for lengh 30 and lengh 50
    * numpy version of logo calculation script


