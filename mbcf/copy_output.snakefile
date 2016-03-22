# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import sys
import os

VIPER_DIR=os.environ.get("VIPER_DIR")
CFCE_RUN=os.environ.get("CFCE_RUN")

out_file_list = []

if not VIPER_DIR:
    print("Execute the snakefile as: VIPER_DIR=/path/to/copy/dir_name snakemake -s copy_output.snakefile")
    sys.exit()

def fusion_output():
    if os.path.isdir("analysis/STAR_Fusion"):
        return VIPER_DIR + "/alignment/gene_fusions/"
    else:
        return ""    


if not CFCE_RUN or int(CFCE_RUN) != 1:
    out_file_list = [
        VIPER_DIR + "/alignment/bam/", 
        VIPER_DIR + "/alignment/bigwig/", 
        fusion_output(), 
        VIPER_DIR + "/diffexp/", 
        VIPER_DIR + "/expression/", 
        VIPER_DIR + "/QC/",  
        VIPER_DIR + "/SNP/", 
        VIPER_DIR + "/plots/"
    ]
else:
    out_file_list = [
        VIPER_DIR + "/analysis/",
        VIPER_DIR + "/viper/"
    ]

rule target:
    input:
        expand("{filename}", filename=out_file_list)


rule copy_cfce_analysis:
    output:
        analysis_dir=VIPER_DIR + "/analysis/"
    shell:
        "cp -rf analysis/* {output.analysis_dir}/"
        " && cp config.yaml " + VIPER_DIR + "/"
        " && cp metasheet.csv " + VIPER_DIR + "/"
        " && cp report.html " + VIPER_DIR + "/$(basename " + VIPER_DIR + ")_report.html" 

rule copy_cfce_viper:
    output:
        viper_dir=VIPER_DIR + "/viper/"
    shell:
        "cp -rf viper/* {output.viper_dir}/"

rule copy_fastqc:
    output:
        fastqc_dir="RNASeq/FastQC/"
    shell: 
        "find fastqc/before_filtering/ -name \"[0-9]*html\" -exec cp -t {output.fastqc_dir} {{}} \;"

rule copy_alignment_bam:
    output:
        bam_dir=VIPER_DIR + "/alignment/bam/"
    shell:
        "find analysis/STAR/ -type f -name \"*\.sorted\.ba*\" -exec cp -t {output.bam_dir} {{}} \;"
        " && cp analysis/STAR/STAR_Align_Report.csv " + VIPER_DIR + "/alignment/Alignment_Report.csv"
        " && cp analysis/STAR/STAR_Align_Report.png " + VIPER_DIR + "/alignment/Alignment_Report.png"


rule copy_alignment_bw:
    output:
        bw_dir=VIPER_DIR + "/alignment/bigwig/"
    shell:
        "find analysis/bam2bw/ -type f -name \"*\.bw\" -exec cp -t {output.bw_dir} {{}} \;"


rule copy_fusion_output:
    output:
        fusion_dir=VIPER_DIR + "/alignment/gene_fusions/"
    shell:
        "for file in $(find analysis/STAR_Fusion/ -type f -name \"*final.abridged\"); do out_file=$(basename $file); tr '\t' ',' <$file 1>{output.fusion_dir}/${{out_file}}.csv; done"

rule copy_diffexp:
    output:
        de_dir=VIPER_DIR + "/diffexp/"
    shell:
        "cp -rf analysis/diffexp/* {output.de_dir}"


rule copy_cufflinks_exp:
    output:
        cuff_dir=VIPER_DIR + "/expression/"
    shell:
        "cp -rf analysis/cufflinks {output.cuff_dir}"
        " && cp analysis/cufflinks/Cuff_Gene_Counts.csv {output.cuff_dir}/Normalized_FPKM_Gene_Counts.csv"
        " && cp analysis/STAR/STAR_Gene_Counts.csv {output.cuff_dir}/Raw_Gene_Counts.csv"


rule copy_RSeQC:
    output:
        qc_dir=VIPER_DIR + "/QC/"
    shell:
        "cp -rf analysis/RSeQC/* {output.qc_dir}"


rule copy_snp:
    output:
        snp_dir=VIPER_DIR + "/SNP/"
    shell:
        "cp -rf analysis/snp/* {output.snp_dir}"


rule copy_summary:
    output:
        plots_dir=VIPER_DIR + "/plots/"
    shell:
        "find analysis/ -type f -name \"*.png\" -exec cp -t {output.plots_dir} {{}} \;"
        " && cp report.html $(basename " + VIPER_DIR + ")_report.html"  
