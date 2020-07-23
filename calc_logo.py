#!/usr/bin/env python
#
# Calculates a sequence logo using Weblogo3 library
# Shows signatures with amino acids on position above certain information threshold (in bits)
#

import string
import random
import sys
import math
from numpy import *

import corebio
#import morecb
from morecb import *
from corebio.seq import Alphabet, Seq, SeqList
from corebio import seq_io
from corebio.utils import isfloat, find_command, ArgumentError
from corebio.moremath import *
from corebio.data import amino_acid_composition
from corebio.seq import unambiguous_rna_alphabet, unambiguous_dna_alphabet, unambiguous_protein_alphabet

from weblogolib import *
import array



### parameters ###
if len(sys.argv)-1 < 3:
    print "Usage:       ./profile_logo.py <matrix filename> <matrix number|all> <signature bit threshold> <cut_from> <cut_to>"
    sys.exit(1)

matrix_filename = sys.argv[1]

if (sys.argv[2] == "all"):
    matrix_number = -1
else:
    matrix_number = int(sys.argv[2])

signature_bit_threshold = float(sys.argv[3])

if len(sys.argv)-1 == 5:
    cut_from = int(sys.argv[4])
    cut_to = int(sys.argv[5])
else:
    cut_from = 1
    cut_to = 50

cut_length = cut_to - cut_from + 1

print "Matrix filename = ", matrix_filename
print "Matrix number = ", matrix_number
print "Signature bit threshold = ", signature_bit_threshold
print "Cut %d positions [from %d to %d]" % (cut_length, cut_from, cut_to)
print

fc = open("composition.csv")
B = zeros((20))
for line in fc:
    (aa, bg) = line.split(",")
    bg = float(bg)/100.0
    B[corebio.seq.unambiguous_protein_alphabet.ords(aa)[0]] = float(bg)
    print aa, float(bg)
fc.close()
print

fl = open("logo_list.html", 'w')
fsignature = open("prototype_signatures.csv", 'w')
fsignature.write("matrix\tsignature\n")

"""
PROTOTYPE 1
BEGIN
SEGMENT VNRKEVSIEGTVTRLWKPSSAAISQVGLIEDESGKTKFTSWVASEQPWIE
MATRIX K=586 N=19305087 P=0.00000050 S=-0.458219 W=2.000000
50     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y
 0 0.092150 0.000000 0.015358 0.039249 0.006826 0.047782 0.001706 0.066553 0.030717 0.018771 0.006826 0.105802 0.027304 0.013652 0.000000 0.148464 0.213311 0.141638 0.001706 0.022184
 1 0.011945 0.003413 0.216724 0.006826 0.010239 0.015358 0.003413 0.025597 0.104096 0.010239 0.001706 0.146758 0.013652 0.006826 0.005119 0.119454 0.247440 0.037543 0.000000 0.013652

"""

freqs = zeros((30, 20))
signature = ""
gap = 0
K = 0
ID=0
format = 0

fin = open(matrix_filename);
for line in fin:
    line = string.strip(line)
    if string.find(line, "MATRIX") != -1:
        (a, b, c) = line.split(" ")
        (t, ID) = b.split("=")
        (t, K) = c.split("=")
        ID = int(ID)
        K = int(K)

        #if (ID < 400):
        #    for x in range(52):
        #       fin.next() # seek forward

        if (matrix_number > 0 and matrix_number != ID):
            for x in range(52):
                fin.next() # seek
        continue
    if string.find(line, "PROTOTYPE") != -1:
        (a, b) = line.split(" ")
        format = int(b)
        continue
    if string.find(line, "SEGMENT") != -1:
        continue
    if string.find(line, "BEGIN") != -1:
        freqs = zeros((50, 20))
        signature = ""
        gap = 0
        continue
    if string.find(line, "END") != -1:

        fsignature.write(str(ID)+"\t"+signature+"\n")

#        data = LogoData.from_counts(corebio.seq.unambiguous_protein_alphabet, freqs, \
#                       prior= ones( (20), float64) /20

        freqs_cut = zeros( (cut_length, 20), float64 )
        for l in range(cut_length):
            for j in unambiguous_protein_alphabet:
                o = unambiguous_protein_alphabet.ords(j)[0]
                freqs_cut[l, o] = freqs[cut_from - 1 + l,o]

        data = LogoData.from_counts(corebio.seq.unambiguous_protein_alphabet, freqs_cut, \
                        prior=B    )

        #data = LogoData.from_counts(corebio.seq.unambiguous_protein_alphabet, freqs, prior=ones((20), float64))


        #ent = zeros(  seq_length, float64)
        #entropy_interval = zeros( (seq_length,2) , float64)
        #return cls(seq_length, alphabet, counts, ent, entropy_interval, weight)

        options = LogoOptions()
        options.logo_title = "Profile " + str(ID) + " (K=" + str(K) + ")"
        options.show_yaxis = True
        options.show_errorbars = False
        options.show_title = True
        #options.logo_start = cut_from
        #options.logo_end = cut_to
        #options.first_index = 0 - (cut_from - 1) + 1

        #options.composition = 'equiprobable'
        options.stacks_per_line = 50
        options.stack_width = std_sizes['large']

        options.creator_text =''

        options.color_scheme = std_color_schemes['chemistry']
        format = LogoFormat(data, options)
        fout = open('logo/' + str(ID) + '_logo.png', 'w')
        #png_formatter( data, format, fout) # png_print_formatter, #eps_formatter
        png_print_formatter( data, format, fout) # png_print_formatter, #eps_formatter
        #    fout_eps = open('pattern/' + str(pattern_id) + '_logo.eps', 'w')
        #    eps_formatter( data, format, fout_eps) # png_print_formatter, #eps_formatter
        fout.close()
        #    fout_eps.close()
        print str(ID), signature #, "   cut:", signature[cut_from-1:cut_to-1]

        continue

    if line.find("50", 0, 2) != -1:
        continue

    pos = int(line[0:2])
    freqline = line.split(" ")
    freqline = freqline[1:]
    aa_index=0
    freqdict = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        freqdict[aa] = freqline[aa_index]
        aa_index += 1

    #print "processing pos", pos, freqdict

    #fl.write("<a href=\"prototype_profiles/" + str(pattern_id) + ".html\"><img src=\"prototype_logo/" + str(pattern_id) + "_logo.png\" width=\"600\"></a><br>\n")

    if (K>0):
        e = (20 - 1) / (2.0 * math.log(2) * K) # e is the small-sample correction for an amino acid alignment of N letters
    else:
        e = 0.0

    k = pos
    H = 0.0
    #print "pos ", pos
    for j in corebio.seq.unambiguous_protein_alphabet:
        o = corebio.seq.unambiguous_protein_alphabet.ords(j)[0]
        #print j, freqdict[j]
        freqs[k, o] = float(string.lstrip(freqdict[j]))
        if (freqs[k, o] > 0):
        #H +=  freqs[k, o] * math.log(freqs[k, o]/B[o])/math.log(2)
            H -=  freqs[k, o] * (math.log(freqs[k, o])/math.log(2))

                                                # H = uncertainty at current position k measured in bits per position
        #R = H # total information at the position is decrease in uncertainty

        R = (math.log(20)/math.log(2)) - H

        #print H, e, R


    variants = []
    #gap = 0
    for j in corebio.seq.unambiguous_protein_alphabet:
        o = corebio.seq.unambiguous_protein_alphabet.ords(j)[0]
        #print 1.0*freqs[k, o]/K * R
        ## calc consensus string taking one bit

        if freqs[k, o] * R >= signature_bit_threshold and j not in variants:
            variants.append(j)

        freqs[k, o] = int(K*freqs[k, o])
        print freqs

    #print K, sum(freqs[k,]), freqs[k,]

    """
    if len(variants) > 0 or k == 49:
        if gap > 0:
            if gap == 1:
                signature += "*"
            else:
                signature += "("+str(gap)+"*)"
            gap = 0
        if len(variants) == 1:
            signature += variants[0]
        elif len(variants) > 0:
            signature += "["+"".join(list(set(variants)))+"]"
    else:
            gap += 1
    """
    #print pos, len(variants)
    if len(variants) > 0 or k == 49:
        if len(variants) == 1:
            signature += variants[0]
        elif len(variants) > 0:
            signature += "["+"".join(list(set(variants)))+"]"
    else:
        signature += "-"

fin.close()
fl.close()
fsignature.close()
