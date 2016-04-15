#!/usr/bin/python
"""Given a UCSC refSeq bed file with refseq IDS for gene names/ids, converts
the refseq IDs to hugo names

Steps for generating refseq.bed:
   a. goto http://genome.ucsc.edu/
   b. click Tables, clade: mammals, genome: Mouse, assembly: mm9,
      groups: Genes and Gene Prediction tracks, track: Refseq Genes
      table: refGene
      output format: bed

Steps for generating link file:
      1. get gene names: dump another table w/ gene names
      do the same as step 1, but in
      output format: selected fields from primary and related tables

      Then:
      get output
      choose name, name2 as the field (along w/ chr, start, end)
      get output

"""
import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -b bed_file -l link_file (that has refseq and hgnc)"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-b", "--bed", help="bed file (to convert)")
    optparser.add_option("-l", "--link", help="link file (contains refseq col0 and hgnc col4)")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.bed or not options.link:
        optparser.print_help()
        sys.exit(-1)

    #read in linking table, dump into dictionary: REF_ID: HUGO NAME
    f = open(options.link)
    dict = {}
    for l in f:
        if l.startswith("#"):
            continue
        else:
            tmp = l.strip().split("\t")
            dict[tmp[0]] = tmp[4].upper()
    f.close()

    #CONVERT:
    f = open(options.bed)
    for l in f:
        if l.startswith("#"):
            continue
        else:
            tmp = l.strip().split("\t")
            refID = tmp[3]
            #REPLACE
            tmp[3] = dict[refID]
            print("\t".join(tmp))
    f.close()

if __name__=='__main__':
    main()

