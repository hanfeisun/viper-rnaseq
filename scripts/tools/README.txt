These are some tools/scripts to help you process files.  This README contains
a brief description about each tool in this dir.

convertRefseq2Hugo_bed.py - Most people will use the UCSC genome table to 
generate the gene bed files.  Sometimes this gene bed files might have 
REFSEQ IDs rather than Hugo IDs/Names and you'd like to convert to using 
HUGO names.  AFTER generating a linking file (a file that has both REFSEQ ID
and HUGO names for each gene--also done through UCSC genome table), you can 
use this script to subsitute in Hugo Names in place of REFSEQ ids.

convertRefseq2Hugo_gtf.py - VERY similar to convertRefseq2Hugo_bed.py BUT for 
gtfs (from UCSC).
