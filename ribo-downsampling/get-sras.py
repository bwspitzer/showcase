#!/usr/bin/env python3

# Script for downloading samples by SRA from the NCBI Gene Expression Omnibus.
# Used in a Snakemake pipeline by 'Snakefile_twopass_v3'.
# See that file for details.

import os
import sys

os.system("mkdir -p sras")
sample=str(sys.argv[1])
prefix = str(sys.argv[2])

os.system("wget -O sras/{2}.sra {0}/{1}/{2}/{2}.sra".format(prefix, sample[:6], sample))
