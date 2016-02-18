#!/usr/bin/env python2.7
"""Script to calculate SNP correlations among samples.

Method Corr: let samples S1 and S2 have a set of common snps C.  M is the 
subset of C where the SNP from S1 = SNP from S2.

Corr = | M | / | C |
"""
import os
import sys
from optparse import OptionParser

def readSNPFile(filename):
    """Given a varscan output file, returns the SNPS on chr6"""
    _readCut = 10
    f = open(filename)
    header = f.readline().strip().split('\t')
    
    snps = {}
    highSnps = {}

    for line in f:
        tmp = line.strip().split("\t")
        #DROP!if tmp[0] == "chr6":
        if tmp[0]:
            #TAKE this SNP!
            snps[tmp[1]] = tmp[3]
            #highly reliable snps are at least 10 reads for each allele and 
            #at least 25% alternative allele
            if ((int(tmp[4]) >= _readCut) and (int(tmp[5]) >= _readCut) and\
                    (float(tmp[6][:-1]) >= 25)):
                highSnps[tmp[1]] = tmp[3]
    f.close()

    return (snps, highSnps)

def main():
    usage = "USAGE: %prog  file1 file2 ... fileN"
    #optparser.OptionParser.format_epilog = lambda self, formatter: self.epilog
    optparser = OptionParser(usage=usage)
    (options, args) = optparser.parse_args(sys.argv)

    if len(args) < 3:
        optparser.print_help()
        sys.exit(-1)

    snp = {}
    snp_hi = {}
    samplesLs = []
    for f in args[1:]:
        #parse out the name: e.g. /some/path/to/SAMPLE.snp.txt
        sample = (f.split("/")[-1]).strip(".snp.txt")
        (snp[sample], snp_hi[sample]) = readSNPFile(f)
        samplesLs.append(sample)

    corrMat = []
    for (i, s1) in enumerate(samplesLs):
        corrMat.append([])
        for (j, s2) in enumerate(samplesLs):
            if i == j:
                corrMat[i].append(1.0)
            else:
                snps1 = snp_hi[s1]
                snps2 = snp[s2]
                common = set(snps1.keys()).intersection(set(snps2.keys()))
                ct = sum([1 for s in common if snps1[s] == snps2[s]])
                #print "%s:%s %.8f" % (s1, s2, float(ct)/(len(common)+.001))
                corr = float(ct)/(len(common)+.001)
                corrMat[i].append(corr)

    #print corrMat
    #header:
    print("\t".join(["Samples"]+samplesLs))
    for (i,row) in enumerate(corrMat):
        print("\t".join([samplesLs[i]]+map(str,row)))
            

if __name__=='__main__':
    main()

