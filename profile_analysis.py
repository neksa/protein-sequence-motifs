#!/usr/bin/env python
#
# Analyses and filters profile matrices
# Shows signatures with amino acids on position above certain information threshold (in bits)
#

import string
import random
import sys
import math
import re
from numpy import *

import array

protein_alphabet = "ACDEFGHIKLMNPQRSTVWY"

### parameters ###
if len(sys.argv)-1 < 2:
    print "Usage:       ./profile_analysis.py <matrix filename> <matrix number|all>"
    sys.exit(1)

matrix_filename = sys.argv[1]

#calc_stats = False
calc_stats = True

matrix_numbers = []

try:
    id_list_file = open(sys.argv[2])
    for line in id_list_file:
        line = line.strip()
        matrix_numbers.append(int(line))
    id_list_file.close()
except IOError:
    print("IO error")
    if (sys.argv[2] != "all"):
        matrix_numbers.append(int(sys.argv[2]))
except ValueError:
    print("Value error")

print "Matrix filename = ", matrix_filename
print "Matrix numbers = ", matrix_numbers

fc = open("composition.csv")
B = zeros((20))
for line in fc:
    (aa, bg) = line.split(",")
    bg = float(bg)/100.0
    B[protein_alphabet.index(aa)] = float(bg)
    print aa, float(bg)
fc.close()
print

freqs = zeros((50, 20))
signature = ""
gap = 0
K = 0
ID=0
format = 0
matrix_counter = 0
skip = False
signature_bit_threshold = 1

fin = open(matrix_filename);
fout = open("filtered.matrix", "w");

if calc_stats:
    fmatrix_stats = open("matrix.stats.tab", "w")
    fmatrix_stats.write("id\tK\tS\tlag\taa\tautocorr\tmin_pos_inf\tmean_pos_inf\tmax_pos_inf\tsignature\n")

max_pos_inf = 0.0
min_pos_inf = 100
mean_pos_inf = 0.0
max_autocorrelation = 0.0
skipped_counter = 0

matrix_buf = ""

for line in fin:
    matrix_buf += line
    line = string.strip(line)

    if string.find(line, "END") != -1:
        #if (min_pos_inf > 0.6):
        #       print "Signature:", str(ID)+"\t"+signature+"\n"
        #        print "Max positional information:", max_pos_inf

        if calc_stats:
            max_lag_r = 0.0;
            max_lag = 0
            max_lag_aa = "A"
            for lag in range(2,16):
                n = 50.0 - lag

                for o in range(len(protein_alphabet)):
                    lag_r_o = 0.0;
                    f_mean = t_mean = 0.0
                    for i in range(50-lag):
                        f_mean += freqs[i,o]
                        t_mean += freqs[i+lag,o]
                    f_mean /= n
                    t_mean /= n

                    f_sigma2 = t_sigma2 = 0.0
                    for i in range(50-lag):
                        f_sigma2 += (freqs[i,o] - f_mean) * (freqs[i,o] - f_mean)
                        t_sigma2 += (freqs[i+lag,o] - t_mean) * (freqs[i+lag,o] - t_mean)
                    f_sigma = math.sqrt(f_sigma2/n)
                    t_sigma = math.sqrt(t_sigma2/n)

                    for i in range(50-lag):
                        lag_r_o += (freqs[i,o] - f_mean) * (freqs[i+lag,o] - t_mean)
                    lag_r_o /= f_sigma * t_sigma * (n-1)
                    lag_r_o = round(lag_r_o, 2)

                    #print "lag=",lag,"R=",lag_r_o
                    if lag_r_o > max_lag_r:
                        max_lag_r = lag_r_o
                        max_lag = lag
                        max_lag_aa = protein_alphabet[o]

            print "%d\t%d\t%.2f\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%s"% (ID, K, S, max_lag, max_lag_aa, max_lag_r, min_pos_inf, mean_pos_inf, max_pos_inf, signature)

            #ID matrix_counter
            fmatrix_stats.write("%d\t%d\t%.2f\t%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n"% (ID, K, S, max_lag, max_lag_aa, max_lag_r, min_pos_inf, mean_pos_inf, max_pos_inf, signature))

            #548    101     -0.32   13      H       0.50    0.56    1.73    3.36

            #if max_lag_r < 0.9 and max_pos_inf >= 4 and mean_pos_inf < 1.6 and min_pos_inf < 0.4 and K < 900:
            if max_lag_r < 0.95 and max_pos_inf >= 3 and mean_pos_inf < 2.5 and min_pos_inf < 1 and K < 1000 and K > 30:
                fout.write(matrix_buf)
            else:
                print "SKIP by statistics filter"
                skipped_counter += 1
        else:
            if not skip:
                fout.write(matrix_buf)
        matrix_counter += 1
        matrix_buf = ""
        skip = False
        continue

    if skip: continue

    if line.find("50", 0, 2) != -1:
        continue

    if string.find(line, "MATRIX") != -1:
        S = 0.0
        K = 0
        ID = 0

        pairs = line.split(" ")
        for pair in pairs:
            if pair.find("=") != -1:
                (var, value) = pair.split("=", 1)
                if var == "ID": ID = int(value)
                if var == "K": K = int(value)
                if var == "S": S = float(value)

        max_pos_inf = 0.0
        min_pos_inf = 100
        mean_pos_inf = 0

        if (TYPE == 1 or TYPE == 3):
            ID = matrix_counter

        if (len(matrix_numbers) and ID not in matrix_numbers):
            skip = True
            #print "SKIPPING..", ID
            #for x in range(52):
            #   matrix_buf = ""
            #   fin.next() # seek
        continue
    if string.find(line, "PROTOTYPE") != -1 or string.find(line, "PROFILE") != -1:
        typeline = re.split("\s+", line)
        TYPE = int(typeline[1])
        skip = False
        if (TYPE not in (1,2,3,4)):
            print "Unsupported matrix type: ", line
            skip = True
        continue
    if string.find(line, "SEGMENT") != -1:
        continue
    if string.find(line, "BEGIN") != -1:
        freqs = zeros((50, 20))
        signature = ""
        gap = 0
        continue

    ######### process next matrix position ##########

    pos = int(line[0:2])
    k = pos

    freqline = re.split("\s+", line)
    freqline = freqline[1:]
    aa_index = 0
    for o in range(len(protein_alphabet)):
        #print k, o, aa_index, freqline
        freqs[k, o] = float(freqline[aa_index].lstrip())
        aa_index += 1
    #print freqs

    if (K>0):
        e = (20 - 1) / (2.0 * math.log(2) * K) # e is the small-sample correction for an amino acid alignment of N letters
    else:
        e = 0.0

    H = 0.0

    for o in range(len(protein_alphabet)):
        if (TYPE==3 or TYPE==1):
            freqs[k, o] /= K

        if (freqs[k, o] > 0):
            #H +=  freqs[k, o] * math.log(freqs[k, o]/B[o])/math.log(2)
            H -=  freqs[k, o] * (math.log(freqs[k, o])/math.log(2))
        # H = uncertainty at current position k measured in bits per position
        #R = H # total information at the position is decrease in uncertainty
        R = (math.log(20)/math.log(2)) - H


    if max_pos_inf < R: max_pos_inf = R
    if min_pos_inf > R: min_pos_inf = R
    mean_pos_inf += R/50.0

    variants = []
    for o in range(len(protein_alphabet)):
        j = protein_alphabet[o]
        #print 1.0*freqs[k, o]/K * R
        ## calc consensus string taking one bit
        if freqs[k, o] * R >= signature_bit_threshold and j not in variants:
            variants.append(j)
        freqs[k, o] = int(K*freqs[k, o])

    if len(variants) > 0 or k == 49:
        if len(variants) == 1:
            signature += variants[0]
        elif len(variants) > 0:
            signature += "["+"".join(list(set(variants)))+"]"
    else:
        signature += "-"

fin.close()
fout.close()
if calc_stats:
    fmatrix_stats.close()

print "DONE: processed", matrix_counter, "profiles,", "removed", skipped_counter, "profiles"
