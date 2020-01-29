#!/usr/bin/env python3

# Script for collecting splice junctions for an aligned set of samples, for the purpose of 
#   performing a second alignment with STAR.
# Used in a Snakemake pipeline by 'Snakefile_twopass_v3'.
# See that file for details.

import pandas as pd
import csv
import sys
import os

os.system("mkdir -p SJ")
INFILE = "first_pass/{0}SJ.out.tab".format(sys.argv[1])
OUTFILE = "SJ/{0}_tSJ.out.tab".format(sys.argv[1])

df = pd.read_csv(INFILE,sep='\t',header=None,names = ['chr',
                                                      'intr_1st',
                                                      'intr_last',
                                                      'strand',
                                                      'motif',
                                                      'annotated',
                                                      'unique_crossing',
                                                      'multi_crossing',
                                                      'max_overhang'])

df = df.loc[(df['chr'] != 'chrM') &
            (df['motif'] > 0) &
            (df['annotated'] > 0) &
            (df['unique_crossing'] > 2)]

df.to_csv(OUTFILE, sep ='\t', header = False, index = False)


