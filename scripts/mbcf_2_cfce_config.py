#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

# @AUTHOR: Mahesh Vangala
# @Email: vangalamaheshh@gmail.com
# @Date: Apr,8,2016

import argparse
import sys
import os.path
import yaml
import pandas

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--metasheet', required=True, help="Provide a mbcf specific metasheet: Will have \
        'Filename\tSampleName\tCondition1\tConditionN\tcomp_1vs2\tcomp_1vsN'")
    parser.add_argument('-r','--reference', required=True, help="Reference Organism. Supported values are [hg19,mm9]")
    args = parser.parse_args()
    return args

def getYaml( org ):
    org_yaml = './viper/mbcf/' + org + '.yaml'
    with open(org_yaml, "r") as fh:
        return yaml.safe_load(fh)

def getSampleInfo( metasheet ):
    df = pandas.read_csv( metasheet,index_col=1 )
    return df[df.columns[0]].to_dict()  

def getSampleFiles( sampleInfo ):
    samples = {}
    for samplename, leftmate in sampleInfo.items():
        rightmate = leftmate.replace("_R1_","_R2_")
        if os.path.isfile("./data/" + rightmate):
            samples[samplename] = ["./data/" + leftmate, "./data/" + rightmate]
        else:
            samples[samplename] = ["./data/" + leftmate]
    
    return samples

def print_config_yaml( configObj ):
    with open( "config.yaml", "w" ) as outfile:
        outfile.write( yaml.dump( configObj, default_flow_style=False ) )

def print_metasheet( metasheet ):
    df = pandas.read_csv( metasheet )
    with open("metasheet.csv","w") as outfile:
        outfile.write(df[df.columns[1:]].to_csv(index=False))

def symlink_ref_yaml( org ):
    os.symlink( './viper/mbcf/' + org + '_ref.yaml', 'ref.yaml' )

args = parseArgs()
if args.reference not in ["hg19","mm9"]:
    sys.stderr.write( args.reference + " is not supported. Please provide one of the values from [hg19,mm9]. Exiting ...!\n" )
    sys.exit(1)

if not os.path.isfile(args.metasheet):
    sys.stderr.write( args.metasheet + " file doesn't exist. Exiting ...!\n" )
    sys.exit(1)

 
if __name__ == '__main__':
    configObj = getYaml( args.reference )
    sampleInfo = getSampleInfo( args.metasheet )  
    sampleFiles = getSampleFiles( sampleInfo )
    configObj["the_samples"] = sampleFiles
    print_config_yaml( configObj )
    print_metasheet( args.metasheet )
    symlink_ref_yaml( args.reference )
