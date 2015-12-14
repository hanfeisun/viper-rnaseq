# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule target:
    input:
        "RNASeq/FastQC/", "RNASeq/Alignment/", "RNASeq/Expression/FPKM/Genes/", "RNASeq/Expression/FPKM/Transcripts/",  "RNASeq/Expression/Raw/", "RNASeq/Summary/"

rule copy_fastqc:
    output:
        fastqc_dir="RNASeq/FastQC/"
    shell: 
        "find fastqc/before_filtering/ -name \"[0-9]*html\" -exec cp -t {output.fastqc_dir} {{}} \;"

rule copy_alignment:
    output:
        align_dir="RNASeq/Alignment/",
        bam_dir="RNASeq/Alignment/bam/",
        big_wig="RNASeq/Alignment/bigwig/"
    shell:
        "find analysis/STAR/ -type f -name \"*\.sorted\.ba*\" -exec cp -t {output.bam_dir} {{}} \;"
        " && cp analysis/STAR/STAR_Gene_Counts.csv {output.align_dir}/Raw_Gene_Counts.csv"
        " && cp analysis/STAR/STAR_Align_Report.csv {output.align_dir}/Alignment_Report.csv"
        " && find analysis/bam2bw/ -type f -name \"*\.bw\" -exec cp -t {output.big_wig} {{}} \;"


rule copy_fpkm_genes:
    output:
        fpkm_gene_dir="RNASeq/Expression/FPKM/Genes/"
    shell:
        "find analysis/cufflinks/ -type f -name \"*\.genes\.fpkm_tracking\" -exec cp -t {output.fpkm_gene_dir} {{}} \;"
        " && cp analysis/cufflinks/Cuff_Gene_Counts.csv RNASeq/Expression/FPKM_Gene_Counts.csv"

rule copy_fpkm_transcripts:
    output:
        fpkm_trans_dir="RNASeq/Expression/FPKM/Transcripts/"
    shell:
        "find analysis/cufflinks/ -type f -name \"*\.isoforms\.fpkm_tracking\" -exec cp -t {output.fpkm_trans_dir} {{}} \;"

rule copy_raw:
    output:
        raw_dir="RNASeq/Expression/Raw/"
    shell:
        "find analysis/STAR/ -type f -name \"*\.counts\.tab\" -exec cp -t {output.raw_dir} {{}} \;"
        " && cp analysis/STAR/STAR_Gene_Counts.csv RNASeq/Expression/Raw_Gene_Counts.csv"

rule copy_summary:
    output:
        "RNASeq/Summary/"
    shell:
        "cp analysis/STAR/STAR_Align_Report.png {output}/Alignment_Report.png"
        " && cp analysis/STAR/STAR_Gene_Counts.csv {output}/Raw_Gene_Counts.csv"
        " && cp analysis/cufflinks/Cuff_Gene_Counts.csv {output}/FPKM_Gene_Counts.csv"
        " && cp analysis/RSeQC/read_distrib/read_distrib.png {output}/"
        " && cp analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png {output}/"
        " && cp analysis/RSeQC/gene_body_cvg/geneBodyCoverage.curves.png {output}/"
        " && find . -mindepth 1 -maxdepth 1 -type f -name \"*csv\" -a ! -name \"\.*\" -exec cp -t {output} {{}} \;"
