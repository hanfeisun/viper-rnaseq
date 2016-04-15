#!/usr/bin/python
"""Given a UCSC refSeq GTF file with refseq IDS for gene ids, converts
the refseq IDs to hugo names

Steps for generating refseq.gtf:
   a. goto http://genome.ucsc.edu/
   b. click Tables, clade: mammals, genome: Mouse, assembly: mm9,
      groups: Genes and Gene Prediction tracks, track: Refseq Genes
      table: refGene
      output format: gtf

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
    usage = "USAGE: %prog -g gtf_file -l link_file (that has refseq and hgnc)"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-g", "--gtf", help="gtf file (to convert)")
    optparser.add_option("-l", "--link", help="link file (contains refseq col0 and hgnc col4)")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.gtf or not options.link:
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
    f = open(options.gtf)
    for l in f:
        if l.startswith("#"):
            continue
        else:
            tmp = l.strip().split("\t")
            #select out last col -- gene attributes list
            attribs = tmp[-1].strip().split(";")
            #find the gene id attrib
            refID = filter(lambda x: x.startswith('gene_id'), attribs)
            if refID:
                refID = refID[0].split(" ")[1]
                #NOTE: we get something like "NM_001195025"; need to eval as
                refID = eval(refID)
                
                #DROP the old gene_id attrib
                attribs = [a for a in attribs if not a.startswith("gene_id")]
                #INSERT the new gene_id at the head of the list
                attribs.insert(0, "gene_id \"%s\"" % (dict[refID]))
                
                #INSERT back into tmp
                tmp[-1] = ";".join(attribs)

            print("\t".join(tmp))
    f.close()

if __name__=='__main__':
    main()

