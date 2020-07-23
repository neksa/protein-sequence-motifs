#!/usr/bin/env python
# Profile to Logo
# Alexander Goncearenco 2009
#
# Calculates a sequence logo using Weblogo3 library
# Shows signatures with amino acids on position above certain information threshold (in bits)
#
# accepts frequency matrices with fixed length (format PROTOTYPE 2)
# will look for composition.csv in current path (the file is in examples directory)
#
# Dependencies: numpy, weblogolib, corebio
#

import string
import random
import sys
import math
import re
import os
from numpy import *

import corebio
#import morecb
#from morecb import *
from corebio.seq import Alphabet, Seq, SeqList
from corebio import seq_io
from corebio.utils import isfloat, find_command, ArgumentError
from corebio.moremath import *
from corebio.data import amino_acid_composition
from corebio.seq import protein_alphabet, unambiguous_protein_alphabet

from weblogolib import *
import array

### parameters ###
print sys.argv

if len(sys.argv)-1 < 3:
    print "Usage: ./profile_logo.py <matrix filename> <matrix number|all|list file name> <signature bit threshold> <cut_from> <cut_to>"
    print "Requires: composition.csv"
    print "Creates profile_signatures.csv profile_list.html and directory ./profile_logo "
    sys.exit(1)

matrix_filename = sys.argv[1]

matrix_number = []

if (sys.argv[2] != "all" and not(sys.argv[2].isdigit())):
    flist = open(sys.argv[2])
    for line in flist:
        matrix_number.append(int(line.strip()))
    flist.close()
else:
    if (sys.argv[2] != "all"):
        matrix_number.append(int(sys.argv[2]))

mode = 2 #1

signature_bit_threshold = float(sys.argv[3])

print_quality = False
if len(sys.argv)-1 >= 4:
    if (sys.argv[4]) == "print":
        print_quality = True

cut_specified = False
if len(sys.argv)-1 >= 6:
    cut_from = int(sys.argv[5])
    cut_to = int(sys.argv[6])
    cut_specified = True
else:
    cut_from = 1
    cut_to = 50

cut_length = cut_to - cut_from + 1

print "Matrix filename = ", matrix_filename
print "Matrix number ", len(matrix_number)
print "Signature bit threshold = ", signature_bit_threshold
print "Length max %d positions [from %d to %d]" % (cut_length, cut_from, cut_to)
print

if not os.path.exists("profile_logo"):
    os.makedirs("profile_logo")

if not os.path.exists("profile_logo/print"):
    os.makedirs("profile_logo/print")

if os.path.exists("composition.csv"):
    fc = open("composition.csv")
    B = zeros(len(protein_alphabet))
    for line in fc:
        (aa, bg) = line.split(",")
        bg = float(bg)/100.0
        B[protein_alphabet.ords(aa)[0]] = float(bg)
        print aa, float(bg)
    fc.close()
    print
else:
    B = zeros(len(protein_alphabet))
    # Attention: Archaeal frequencies by default
    B[protein_alphabet.ords('A')[0]] = 0.7846065
    B[protein_alphabet.ords('R')[0]] = 0.5395213
    B[protein_alphabet.ords('N')[0]] = 0.3854661
    B[protein_alphabet.ords('D')[0]] = 0.5671408
    B[protein_alphabet.ords('C')[0]] = 0.0983962
    B[protein_alphabet.ords('E')[0]] = 0.7457194
    B[protein_alphabet.ords('Q')[0]] = 0.2283004
    B[protein_alphabet.ords('G')[0]] = 0.7429558
    B[protein_alphabet.ords('H')[0]] = 0.1708365
    B[protein_alphabet.ords('I')[0]] = 0.7471997
    B[protein_alphabet.ords('L')[0]] = 0.9529720
    B[protein_alphabet.ords('K')[0]] = 0.5845627
    B[protein_alphabet.ords('M')[0]] = 0.2372575
    B[protein_alphabet.ords('F')[0]] = 0.3902878
    B[protein_alphabet.ords('P')[0]] = 0.4283092
    B[protein_alphabet.ords('S')[0]] = 0.6101052
    B[protein_alphabet.ords('T')[0]] = 0.5260790
    B[protein_alphabet.ords('W')[0]] = 0.1027624
    B[protein_alphabet.ords('Y')[0]] = 0.3727149
    B[protein_alphabet.ords('V')[0]] = 0.7847484

#fl = open("profile_list.html", 'w')
fsignature = open("profile_signatures.csv", 'w')
fsignature.write("matrix\tsignature\n")

"""
PROTOTYPE 1
BEGIN
SEGMENT HVHPKDLISGEITPIERRGYPAPIVNHNLRQKQFKALYNQLKAAIAEPEA
MATRIX K=637 N=19305087 P=0.00000005 S=0.708697 W=1.000000
50     A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y
 0    57    11     2     6     9     4     4    91    13    20    11     1     3     5     4    54   166   161     1    14
 1    37     9    15    17     1     9    13    92    36    38     6    28    23     7    36    52    35   149     7    27
 2    43     4    65   129     2    20    10    13    22    12     3    25   111    11    18    78    46    16     0     9

PROTOTYPE 2
BEGIN
MATRIX ID=0 K=637
50        A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
 0 0.089482 0.017268 0.003140 0.009419 0.014129 0.006279 0.006279 0.142857 0.020408 0.031397 0.017268 0.001570 0.004710 0.007849 0.006279 0.084772 0.260597 0.252747 0.001570 0.021978
 1 0.058085 0.014129 0.023548 0.026688 0.001570 0.014129 0.020408 0.144427 0.056515 0.059655 0.009419 0.043956 0.036107 0.010989 0.056515 0.081633 0.054945 0.233909 0.010989 0.042386

PROFILE 3
BEGIN
MATRIX K=45 L=2
   2    A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    .
   0    1    0    0    2    0    1    1    1    1    0    1    2    0    0    0    0    0    1    0    1   33
   1    1    0    0    0    2    0    0    0    1    1    0    1    2    0    1    0    2    1    0    0   33

PROFILE 4
BEGIN
MATRIX ID=0 K=30 L=30
30        A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
 0 0.089482 0.017268 0.003140 0.009419 0.014129 0.006279 0.006279 0.142857 0.020408 0.031397 0.017268 0.001570 0.004710 0.007849 0.006279 0.084772 0.260597 0.252747 0.001570 0.021978
 1 0.058085 0.014129 0.023548 0.026688 0.001570 0.014129 0.020408 0.144427 0.056515 0.059655 0.009419 0.043956 0.036107 0.010989 0.056515 0.081633 0.054945 0.233909 0.010989 0.042386

"""

freqs = zeros((30, len(protein_alphabet))) # just to init
signature = ""
gap = 0
K = 0
L = 0
ID = 0
TYPE = 0
skip = True
matrix_counter = 0
first_meaningful_position = -1
last_meaningful_position = -1
origin = ""

fin = open(matrix_filename);
for line in fin:
    if string.find(line, "PROTOTYPE") != -1 or string.find(line, "PROFILE") != -1:
        matrix_counter += 1
        typeline = re.split("\s+", line)
        TYPE = int(typeline[1])
        skip = False
        if (TYPE not in (1,2,3,4)):
            print "Unsupported matrix type: ", line
            skip = True
        continue

    if string.find(line, "SEGMENT") != -1:
        continue

    if string.find(line, "ORIGIN") != -1:
        if len(line) <= 8:
            continue
        (_, origin) = re.split("\s+", line.strip(), 1)
        continue

    if string.find(line, "BEGIN") != -1:
        #freqs = zeros((50, 20))
        signature = ""
        gap = 0
        first_meaningful_position = -1
        last_meaningful_position = -1
        continue

    line = string.strip(line)
    if string.find(line, "MATRIX") != -1:
        matrixline = re.split("\s+", line)
        for x in matrixline[1:]:
            (t, v) = x.split("=")
            if (t == "ID"):
                ID = int(v)
            if (t == "K"):
                K = int(v)
            if (t == "L"):
                L = int(v)

        if (TYPE < 3):
            L = 50

        cut_length = min(cut_to, L) - min(cut_from, L) + 1

        if (TYPE == 1 or TYPE == 3):
            ID = matrix_counter - 1

        print 'TYPE=%d, ID=%d, K=%d, L=%d [%d:%d-%d]'% (TYPE, ID, K, L, cut_length, min(cut_from, L), min(cut_to, L))
        freqs = zeros((L, len(protein_alphabet)));

        if (len(matrix_number) and (ID not in matrix_number)):
            # skip matrix length + 1 header line + 1 end line
            for x in range(L+2): fin.next()
        continue

    if string.find(line, "END") == 0:

        #print freqs

        if (skip): continue

        fsignature.write(str(ID)+"\t"+signature+"\n")

        freqs_cut = zeros( (cut_length, len(protein_alphabet)), float64 )
        B = zeros(len(protein_alphabet))

        for l in range(cut_length):
            for j in protein_alphabet:
                o = protein_alphabet.ords(j)[0]
                freqs_cut[l, o] = freqs[cut_from - 1 + l,o]

        # load frequencies into weblogo compatible object

        P = zeros(len(protein_alphabet))
        #P_freqs = ze
        #print 'len p', len(unambiguous_protein_alphabet), freqs_cut.shape, B.shape
        for j in protein_alphabet:
            B_j = protein_alphabet.ords(j)[0]
            P_j = protein_alphabet.ords(j)[0]
            P[P_j] = B[B_j]

        #data = LogoData.from_counts(protein_alphabet, freqs_cut, prior=B)
        data = LogoData.from_counts(corebio.seq.protein_alphabet, freqs_cut, prior=ones((len(protein_alphabet)), float64)/len(protein_alphabet))

        options = LogoOptions()
        options.logo_title = "Profile " + str(ID) + " (K=" + str(K) + ")" + " " + origin
        options.show_yaxis = True
        options.show_errorbars = False
        options.show_title = True
        #options.logo_start = 10
        #options.logo_end = 40

        #options.composition = 'equiprobable'
        options.stacks_per_line = 50
        options.stack_width = std_sizes['large']
        options.show_fineprint = False
        options.logo_margin = 0.5
        options.creator_text =''

        options.color_scheme = std_color_schemes['chemistry']
        format = LogoFormat(data, options)

        fout = open('profile_logo/' + str(ID) + '_logo.png', 'w')
        png_formatter( data, format, fout)
        fout.close()

        if print_quality:
            options_print = LogoOptions()
            options_print.show_yaxis = True
            options_print.show_xaxis = True
            options_print.show_errorbars = False
            options_print.show_title = False

            if not cut_specified:
                first_index = first_meaningful_position - 1
                #print "CUT: ",first_meaningful_position, last_meaningful_position
            #print "CUT2: ", (0 - (first_index - 1)), (first_meaningful_position - first_index), (last_meaningful_position - first_index)

                options_print.first_index = (0 - (first_index - 1))
                options_print.logo_start = (first_meaningful_position - first_index)
                options_print.logo_end = (last_meaningful_position - first_index)

            options_print.stacks_per_line = 50
            options_print.stack_width = std_sizes['large']
            options_print.creator_text =''
            options_print.color_scheme = std_color_schemes['chemistry']
            options_print.show_fineprint = False
            options_print.logo_margin = 0.5
            format_print = LogoFormat(data, options_print)

            fout = open('profile_logo/print/' + str(ID) + '_logo.png', 'w')
            png_print_formatter( data, format_print, fout)
            fout.close()
            fout = open('profile_logo/print/' + str(ID) + '_logo.eps', 'w')
            eps_formatter( data, format_print, fout)
            fout.close()

        print str(ID), signature
        continue

    if re.match("^%d" % L, line):
        continue

    freqline = re.split("\s+",line)
    #print(freqline)

    pos = int(freqline[0])
    freqline = freqline[1:]
    aa_index=0
    freqdict = {}

    matrix_alphabet = "ACDEFGHIKLMNPQRSTVWY"

    #if (TYPE == 3):
    #   matrix_alphabet += '.'

    for aa in matrix_alphabet:
        #if (aa == '.'):        aa = "-"
        freqdict[aa] = freqline[aa_index]
        aa_index += 1

    #print freqdict

    #print "processing pos", pos, freqdict
    #fl.write("<a href=\"prototype_profiles/" + str(pattern_id) + ".html\"><img src=\"prototype_logo/" + str(pattern_id) + "_logo.png\" width=\"600\"></a><br>\n")

    k = pos
    H = 0.0
    R = 0.0
    #print "pos ", pos
    for j in protein_alphabet:
        o = protein_alphabet.ords(j)[0]

        #print j, o

        if not freqdict.has_key(j): continue

        freqs[k, o] = float(string.lstrip(freqdict[j]))

        if (TYPE==3 or TYPE==1):
            freqs[k, o] /= K

        if (freqs[k, o] > 0):
            #H +=  freqs[k, o] * math.log(freqs[k, o]/B[o])/math.log(2)
            H -=  freqs[k, o] * (math.log(freqs[k, o])/math.log(2))

        #H = uncertainty at current position k measured in bits per position
        R = (math.log(20)/math.log(2)) - H

    # R > 1.5:
    if R > signature_bit_threshold:
        #print "R > 2", k, R
        if first_meaningful_position == -1:
            first_meaningful_position = pos+1
        last_meaningful_position = pos+1

    variants = []
    #gap = 0
    for j in protein_alphabet:
        o = protein_alphabet.ords(j)[0]

        if freqs[k, o] * R >= signature_bit_threshold and j not in variants:
            variants.append(j)

        freqs[k, o] = int(K*freqs[k, o])

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
    if len(variants) > 0 or k == L-1:
        if len(variants) == 1:
            signature += variants[0]
        elif len(variants) > 0:
            signature += "["+"".join(list(set(variants)))+"]"
    else:
        signature += "-"

fin.close()
#fl.close()
fsignature.close()
