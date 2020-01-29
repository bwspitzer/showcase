#!/usr/bin/env python
# coding: utf-8

# This script takes all high-quality transcripts for an annotated genome and
#   divides genes into regions of three types: 1) those which all high-quality
#   transcripts consider to be introns; 2) those which all high-quality
#   transcripts consider to be exons; 3) those which some, but not all, such
#   transcripts consider to be exons.
#
# The establishment of these regions is required for a subsequent algorithm to 
#   estimate how many reads in the latter two types of regions would likely have
#   been retained had the library been prepared through poly-A enrichment rather
#   than ribo-depletion. That algorithm is in the script 'query_names_pt2.py'.

import gffutils
import pysam
import pybedtools
from gffutils import helpers
from collections import defaultdict
import intervaltree
import os
from operator import itemgetter
import pickle
import sys

def to_tree (tree, bedtool):
    for ival in bedtool:
        tree[ival.chrom].add(intervaltree.Interval(ival.start, ival.end))
    return tree
  
def int_list():
    return [0,0]

DBFN = sys.argv[1]
db = gffutils.FeatureDB(DBFN)
chrom = sys.argv[2]
max_TSL = 2
intronic = defaultdict(intervaltree.IntervalTree)
exonic = defaultdict(intervaltree.IntervalTree)
sometimes = defaultdict(intervaltree.IntervalTree)

genes = []
for j in db.region(region=chrom, completely_within = True):
    if j.featuretype == 'gene' and j.attributes['gene_type'][0] == 'protein_coding':
          # gets each 'gene' feature in the specified region
        gene = j.id
        txpts = list(db.children(gene, featuretype='transcript'))    # a list of all transcripts in the 
                                                                     #   annotation for that gene
        temp_txpts = []
        for t in txpts:
            try:    # if this transcript has an entry for 'transcript_support_level'
                    #   and if the level is below the acceptable threshold:
                    #   keep that transcript
                if int([i[1] for i in t.attributes.items() 
                                          if i[0]=='transcript_support_level'][0][0]) <= max_TSL:
                    temp_txpts.append(t)
            except:
                pass

        txpts = temp_txpts
        if (len(txpts) > 0):    # if there are any transcripts with a sufficiently low TSL:
            all_exons = (pybedtools.BedTool([helpers.asinterval(i) for i in db.children(gene, featuretype='exon')]))
            all_exons = all_exons.sort().merge()    # define the ends of the genic region by using the first
                                                    #   and last exon in the annotation as the limits
            gene_extent = pybedtools.BedTool([pybedtools.cbedtools.Interval(chrom = chrom,start = 
                                                                            min(i.start for i in all_exons),end = 
                                                                            max(i.end for i in all_exons))])
            t_introns = []
            for t in txpts:    # for each transcript, get all of the exons
                t_exons = (pybedtools.BedTool([helpers.asinterval(i) for i in db.children(t, featuretype='exon')]))
                t_extent = pybedtools.BedTool([pybedtools.cbedtools.Interval(chrom = chrom,start = min(i.start for i in t_exons), 
                                    end = max(i.end for i in t_exons))])
                t_introns.append(t_extent.subtract(t_exons)) # the transcript may not cover the entire genic region;
                                                             #   to get the introns within a specific transcript, subtract
                                                             #   the exons in that transcript from the extent of that transcript
            if len(t_introns) > 1:    # put together the intronic regions from multiple transcripts if necessary
                always_introns = t_introns[0].intersect(t_introns[1:]).sort().merge()    # get regions which are
                                                                                         #   intronic in all transcripts
                ever_introns = t_introns[0].sort().merge()    # get regions which are intronic in any transcript
                for t in t_introns[1:]:
                    ever_introns = ever_introns.cat(t)
            else:    # as above if there's only one transcript to be considered
                always_introns = t_introns[0]
                ever_introns = always_introns
            always_exons = gene_extent.subtract(ever_introns)
            sometimes_exons = ever_introns.subtract(always_introns)    # get exonic and sometimes regions based on set operations
            temp = [(i.start, i.end, 'exon') for i in always_exons]
            temp.extend([(i.start, i.end, 'intron') for i in always_introns])
            temp.extend([(i.start, i.end, 'sometimes') for i in sometimes_exons])
            temp.sort(key=itemgetter(0))
            regions = [j.id,(0, temp[0][0], 'intergenic')]
            regions.extend(temp)
            regions.extend([(temp[-1][1], float('inf'), 'intergenic')])
            genes.append(regions)

pickle_file = 'pickles/' + str(chrom) + '_pickle'
pck = open(pickle_file, 'wb')
pickle.dump(genes, pck)
pck.close()
