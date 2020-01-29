#!/usr/bin/env python
# coding: utf-8

# This algorithm estimates r, the fraction of reads in an exonic region that would be retained had the library
#   been prepared by poly-A-enrichment rather than by ribo-depletion. For an exon flanked by introns, it does
#   so based on the numbers of spliced and unspliced reads at the splice sites on either side of that exon. For
#   cases where there may be alternative splicing, the algorithm is more complex but the principle is the same.
#
# Genes are considered to consist of a sequence of regions: intronic, exonic, and "sometimes" exonic (presumably
#   due to alternative splicing). Splice junctions separate these regions.
#
# A weakness of this approach is that the algorithm cannot estimate r for genes without introns. In this version
#   of this script, a global value for r is calculated for the entire chromosome and applied to those exons for
#   which there are no neighboring splice junctions.

import pysam
from collections import defaultdict
import sys
import pickle
from random import random as rnd
import numpy as np

def calc_r_recursive (d1f, d0f, region, ds, ds_ctr):
    # This function starts at a given intron and downsamples all regions between that intron and the following
    #   one. Recursion is only required when there is possible alternative splicing and multiple coding regions, 
    #   potentially expressed differently than one another, fall between two introns.
    adj = len(ss_rds[region][spliced_out_before]) - len(ss_rds[region][spliced_out_after])
    # starting at the splice site following region [region]:
    # each read with the current region spliced out means that more reads should
    #   be sampled from the next region; each read with the next region spliced out
    #   means that fewer should be sampled from within the next region
    if regions[region+1][r_type] == 'intron':
        # If an intron has been reached, return in the opposite direction (coords desc)
        d1r = max(0, -adj)
        d0r = len(ss_rds[region][no_splice])
    else:
        tot = sum([len(ss_rds[region][i]) for i in range(3)])
        d1r, d0r, ds, ds_ctr = calc_r_recursive (max(0, d1f + adj), 
                                         max(0, tot - (d1f + adj)), 
                                         region + 1, ds, ds_ctr)
            
    if regions[region][r_type] == 'intron':  # back to the intron where recursion started
        return ds, ds_ctr
    # now that r can be estimated in the 5' -> 3' direction, downsample the reads that fall entirely 
    #   within this region
    try:
        r = (d1f + d1r)/(d1f + d0f + d1r + d0r)
        for q_name in within_region[region]:
            if rnd() < r:
                ds.add(read.query_name)
                ds_ctr['kept'] += 2
            else:
                ds_ctr['discarded'] += 2
    except ZeroDivisionError:
        pass    
             
    # now need to do the reverse-direction operation on the preceding splice site;
    # then adjust d1r and d0r and return them
    if regions[region - 1][r_type] == 'intron':
        # if almost back to the start of the recursion, the preceding splice site has already been downsampled
        return 0, 0, ds, ds_ctr
        
    dsd = len(ss_rds[region-1][spliced_out_before]) + len(ss_rds[region-1][spliced_out_after])
    # this many reads crossing this splice site have already been downsampled
    ds_addl = ((d1f + d1r)/2) - dsd
    # because both neighboring regions are sometimes expressed, some number of unspliced
    #   reads crossing the splice site should be retained
    if ds_addl > 0:
        r = ds_addl/(((d0f + d0r)/2) + ds_addl)
        for q_name in ss_rds[region-1][no_splice]:
            if rnd() < r:
                ds.add(read.query_name)
                ds_ctr['kept'] += 2
            else:
                ds_ctr['discarded'] += 2

    adj = len(ss_rds[region-1][spliced_out_before]) - len(ss_rds[region-1][spliced_out_after])
    tot = sum([len(ss_rds[region-1][i]) for i in range(3)])
    return max(0, d1r - adj), max(0, tot-(d1r - adj)), ds, ds_ctr


def end_exons(ds, ds_ctr, region, direction):
    # This function downsamples exons at the ends of a gene, where there is only one
    #   neighboring splice site. Note that, for the first exon in a gene, direction is 
    #   'reverse'; for the last, it's 'forward'
    if direction == 'reverse':
        direction = -1
    else: 
        direction = 0
    adj = (len(ss_rds[region + direction][spliced_out_before]) - 
           len(ss_rds[region + direction][spliced_out_after]))
    d1 = max(0, adj*(1+2*direction)) # reverses sign if direction == -1
    d0 = len(ss_rds[region][no_splice])
    try:
        r = d1/(d1 + d0)
        dir_temp = (2*direction)+1
        for q_name in within_region[region + dir_temp]:
            if rnd() < r:
                ds.add(read.query_name)
                ds_ctr['kept'] += 2
            else:
                ds_ctr['discarded'] += 2
    except ZeroDivisionError:
        pass
    return ds, ds_ctr
    


in_bam = sys.argv[1]
pickle_file = sys.argv[2]
chrom = sys.argv[3]
outfile = sys.argv[4]

#   make list of regions (incl. 3' intergenic and 5' intergenic)    
temp = open(pickle_file, 'rb')
genes = pickle.load(temp)
temp.close()

r_start, r_end, r_type = 0,1,2  # readable labels for the parts of elements in the list of regions

# note that the splice site with index i is at the 3' end of region i (between region i and region i+1)

last_end = 0
with pysam.AlignmentFile(in_bam) as samfile:
    ds_bam = pysam.AlignmentFile(outfile, 'wb', template = samfile)
    ds_later = set()  # this will be downsampled according to the global r for the chromosome
    ds_ctr = {'kept':0, 'discarded':0}
    by_query_name = defaultdict(set)  # collect all of the reads with a given query_name
    for regions in genes:
        regions = regions[1:]
        introns = [i for i in range(len(regions)) if regions[i][r_type] == 'intron']
        if len(introns) > 0:
            ss = [i[1] for i in regions[:-1]]
            ss_rds = dict() # this will hold the query_names of pairs of reads that cross at least one splice site
                            # example: ss_rds[3][spliced_out_before] contains the query_names of paired reads
                            #   which end in region 4 (after splice site 3) but where region 3 (before splice site 3)
                            #   has been spliced out
            for i in range(len(ss)):
                ss_rds[i] = [set(), set(), set()]
                spliced_out_before, spliced_out_after, no_splice = 0,1,2  
                    # readable labels for the sets within ss_rds
            within_region = dict()
            for i in range(len(regions)):
                if regions[i][r_type] in ['exon', 'sometimes']:
                    within_region[i] = set()
            gene_start = regions[0][r_end]
            gene_end = regions[-1][r_start]
            ds = set() # this will be the set of query_names retained after downsampling
            
            # Some reads fall outside of genic regions or are unpaired; these can't be downsampled by the algorithm
            if gene_start > last_end:
                for read in samfile.fetch(chrom, last_end, gene_start):
                    ds_later.add(read.query_name)
                    by_query_name[read.query_name].add(read)
            for read in samfile.fetch(chrom, gene_start, gene_end):
                by_query_name[read.query_name].add(read)
                if not read.is_proper_pair:
                    ds_later.add(read.query_name)
                    
                else:
                    if read.reference_start < read.next_reference_start:  # read occurs before its mate
                        start_coord = read.reference_start
                        end_coord = start_coord + read.template_length
                    else:
                        start_coord = read.next_reference_start  # read occurs after its mate
                        end_coord = read.reference_end
                    start_region = np.searchsorted(ss, start_coord) # the start of the template falls in this region
                    end_region = np.searchsorted(ss, end_coord) # the end of the template falls in this region
                    if (regions[start_region][r_type] != 'intergenic' and 
                        regions[end_region][r_type] != 'intergenic' and
                        'exon' not in [regions[i][r_type] for i in range(start_region+1, end_region)]):
                        # pairs that start or end outside the gene, and pairs where an exon appears to 
                        #   be spliced out, are considered to be mapping errors
                        if start_region == end_region:
                            try:  # the key only exists for exons and sometimes regions
                                within_region[start_region].add(read.query_name)
                            except KeyError:
                                if regions[start_region][r_type] == 'intergenic':
                                    ds_later.add(read.query_name)
                                    # ignore reads that are entirely within introns
                        elif end_region == start_region +1:  # an unspliced read crossing a splice site
                            ss_rds[start_region][no_splice].add(read.query_name)
                            ds_ctr['discarded'] += 2
                        else:
                            try:
                                ds.add(read.query_name)
                                ds_ctr['kept'] += 2
                                ss_rds[start_region][spliced_out_after].add(read.query_name)  
                                    # a spliced read; the region following start_region has been spliced out
                                ss_rds[end_region-1][spliced_out_before].add(read.query_name)
                                    # a spliced read; the region before end_region has been spliced out
                            except:
                                print(start_coord, end_coord, start_region, end_region)
                                sys.exit(0)
                    else:
                        ds_later.add(read.query_name)

            if regions[introns[0]-1][r_type] in ['exon', 'sometimes']:
                temp1, ds_ctr = end_exons(ds, ds_ctr,
                                                 introns[0], direction = 'reverse')
                ds = ds.union(temp1)
            for intr in introns[:-1]: # exclude the final exon for now b/c it's followed by an intergenic region
                temp1, ds_ctr = calc_r_recursive (0, 0, intr, ds, ds_ctr)
                ds = ds.union(temp1)
            if regions[introns[-1]+1][r_type] in ['exon', 'sometimes']:
                temp1, ds_ctr = end_exons(ds, ds_ctr, introns[-1], direction = 'forward')
                ds = ds.union(temp1)
            for q_name in ds:
                for read in by_query_name[q_name]:
                    ds_bam.write(read)
                del by_query_name[q_name]
                    
            if gene_end > last_end:
                last_end = gene_end
    
    try:
        chrom_r = ds_ctr['kept']/(ds_ctr['kept']+ds_ctr['discarded'])
        for q_name in ds_later:
            if rnd() < chrom_r:
                print ('.',end = '')
                for read in by_query_name[q_name]:
                    ds_bam.write(read)
    except:
        pass
        
    ds_bam.close()

