# Protein sequence motif derivation with extreme sensitivity

Alexander Goncearenco (developed in 2007-2012 at the University of Bergen, Norway)

Designed to derive sequence motifs of elementary functional loops (EFL) in proteins.
Supports MPI parallelization.

The package contains:

    * converge - find conserved motifs in sequence
    * clustering - reduce the number of derived motifs by combining the similar ones
    * search - find matches of a predefined motif in sequece database

Utilities:
    * findPermutatons.R - find positional permutations that destroy relative positional distances and hence the natural motif
    * convertMSA.pl - convert multiple sequence alignment to initial motif
    * composition.csv - average prokaryotic proteomic composition
    * profile_logo.py - generate profile logo graphics

Check other branches as well for modifications.

 