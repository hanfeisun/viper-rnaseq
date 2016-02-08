# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
from collections import defaultdict
import pandas as pd

configfile: "config.yaml"
strand_command=""
cuff_command=""
rRNA_strand_command=""

if( config["stranded"] ):
    strand_command="--outFilterIntronMotifs RemoveNoncanonical"
    cuff_command="--library-type " + config["library_type"]
    rRNA_strand_command="--outFilterIntronMotifs RemoveNoncanonical"
else:
    strand_command="--outSAMstrandField intronMotif"
    rRNA_strand_command="--outSAMstrandField intronMotif"


file_info = defaultdict(list)
ordered_sample_list = []
run_fusion = False

with open( config["metasheet"], "r" ) as meta_fh:
    next(meta_fh)
    for line in meta_fh:
        info = line.strip().split(",")
        file_info[info[0]] = config["samples"][info[0]]
        if( len(file_info[info[0]]) == 2 ):
            run_fusion = True
        if info[0] not in ordered_sample_list:
            ordered_sample_list.append(info[0])

if( run_fusion ):
    if( config["stranded"] ):
        strand_command = " --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000"
    else:
        strand_command = " --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif"


def get_fastq(wildcards):
    return file_info[wildcards.sample]


def fusion_output(wildcards):
    fusion_out_files = []
    if run_fusion:
        for sample in file_info.keys():
            fusion_out_files.append( "analysis/STAR_Fusion/" + sample + "/" + sample + ".fusion_candidates.final" )
    return fusion_out_files

def insert_size_output(wildcards):
    insert_size_out_files = []
    if run_fusion:
        for sample in file_info.keys():
            insert_size_out_files.append( "analysis/RSeQC/insert_size/" + sample + "/" + sample + ".histogram.pdf" )
    return insert_size_out_files

def rRNA_metrics(wildcards):
    if config["star_rRNA_index"] is not None:
        return "analysis/STAR_rRNA/STAR_rRNA_Align_Report.csv"

#LEN: read in comparisons
metadata = pd.read_table(config['metasheet'], index_col=0, sep=',')
comparisons = comparison=[c[5:] for c in metadata.columns if c.startswith("comp_")]

rule target:
    input:
        expand( "analysis/cufflinks/{K}/{K}.genes.fpkm_tracking", K=ordered_sample_list ),
        "analysis/STAR/STAR_Align_Report.csv",
        "analysis/STAR/STAR_Align_Report.png",
        "analysis/STAR/STAR_Gene_Counts.csv",
        "analysis/cufflinks/Cuff_Gene_Counts.csv",
        "analysis/plots/pca_plot.pdf",
        "analysis/plots/heatmapSS_plot.pdf",
        "analysis/plots/heatmapSF_plot.pdf",
        expand( "analysis/RSeQC/read_distrib/{sample}.txt", sample=ordered_sample_list ),
        "analysis/RSeQC/read_distrib/read_distrib.png",
        expand( "analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png", sample=ordered_sample_list ),
        "analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        expand( "analysis/RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf", sample=ordered_sample_list ),
#        expand( "analysis/bam2bw/{sample}/{sample}.bw", sample=ordered_sample_list ),
        expand( "analysis/gfold/{sample}/{sample}.read_cnt.txt", sample=ordered_sample_list ),
        expand("analysis/diffexp/{comparison}/{comparison}.deseq.txt", comparison=comparisons),
        expand("analysis/diffexp/{comparison}/{comparison}_volcano.pdf", comparison=comparisons),

        fusion_output,
        insert_size_output,
        rRNA_metrics

rule report:
    input:
        Unique_Reads="analysis/STAR/STAR_Align_Report.png",
        rRNA_Metrics="analysis/STAR_rRNA/STAR_rRNA_Align_Report.png",
        Read_Distribution="analysis/RSeQC/read_distrib/read_distrib.png",
        Genebody_Coverage_Heatmap="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        Genebody_Coverage_Curves="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.curves.png",
        pca_plot_1="analysis/plots/images/pca_plot_1.png",
        pca_plot_2="analysis/plots/images/pca_plot_2.png",
        pca_plot_3="analysis/plots/images/pca_plot_3.png",
        heatmapSF_plot="analysis/plots/images/heatmapSF_plot.png",
        heatmapSS_plot="analysis/plots/images/heatmapSS_plot.png",
        heatmapSS_cluster="analysis/plots/images/heatmapSS_cluster.png",
        volcano_1="analysis/plots/images/AlvMacsTreatvsVeh_volcano.png",
        volcano_2="analysis/plots/images/EpCamTreatvsVeh_volcano.png",
        volcano_3="analysis/plots/images/TRegsTreatvsVeh_volcano.png",
        volcano_4="analysis/plots/images/CD8treatvsVeh_volcano.png",
        volcano_5="analysis/plots/images/NeutrophilsTreatvsVeh_volcano.png"
    output:
        "report.html"
    run:
        from snakemake.utils import report
        from snakemake.report import data_uri
        unique_reads_en = data_uri( input.Unique_Reads )
        rRNA_metrics_en = data_uri( input.rRNA_Metrics )
        read_distrib_en = data_uri( input.Read_Distribution )
        genebody_hm_en = data_uri( input.Genebody_Coverage_Heatmap )
        genebody_cv_en = data_uri( input.Genebody_Coverage_Curves )
        pca_1_en = data_uri( input.pca_plot_1 )
        pca_2_en = data_uri(input.pca_plot_2)
        pca_3_en = data_uri(input.pca_plot_3)
        heatmapSF_en = data_uri( input.heatmapSF_plot )
        heatmapSS_en = data_uri( input.heatmapSS_plot )
        heatmapSSC_en = data_uri(input.heatmapSS_cluster)
        volcano_1_en = data_uri(input.volcano_1)
        volcano_2_en = data_uri(input.volcano_2)
        volcano_3_en = data_uri(input.volcano_3)
        volcano_4_en = data_uri(input.volcano_4)
        volcano_5_en = data_uri(input.volcano_5)
        report("""
==============================================
VIPER: Visualization Pipeline for RNAseq
==============================================


Alignment Summary
=================
    Samples are aligned against reference organism, using `STAR software`_.

    .. _STAR software: https://github.com/alexdobin/STAR


    The **uniquely mapped read counts** and the **total read counts** for all the samples are summarized in the following image.

    .. image:: {unique_reads_en}

Library Prep Quality Metrics
=============================
Read Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^
    Percentage of reads/sample mapping to specific genomic features

    .. image:: {read_distrib_en}

rRNA removal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Percentage of reads mapping to rRNA sequences. For mRNAseq prep methods ~5% rRNA read mapping is typical. It should be <10% for totalRNAseq. Higher percentages may indicate poor mRNA enrichment or rRNA depletion.

    .. image:: {rRNA_metrics_en}


Genebody Coverage
^^^^^^^^^^^^^^^^^
    For accurate gene expression quantification, mapped reads should be evenly distributed across genebodies.
    Significantly skewed profiles (5' or 3') may introduce quantification bias and/or represent poor quality library preparation.\n

    **Line Plot**

    .. image:: {genebody_cv_en}

    **Heatmap**\n
    This will help you identify biased samples \(lots of blue\)

    .. image:: {genebody_hm_en}


Experimental Design QC
=======================

Principle Component Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Human readable brief description of PCA goes here `wikipedia`_.

    .. _wikipedia: https://en.wikipedia.org/wiki/Principal_component_analysis

    .. image:: {pca_1_en}

    .. image:: {pca_2_en}

    What am I missing if I only look at principle component one and two?

    .. image:: {pca_3_en}


Sample-to-Sample Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. image:: {heatmapSS_en}



Sample-Feature Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. image:: {heatmapSF_en}

    What are *these* genes? \(*link to table coming soon*\)

Differential Gene expression
============================

Volcanos!
^^^^^^^^^
    Brief description goes here

    .. image:: {volcano_1_en}

    .. image:: {volcano_2_en}

    .. image:: {volcano_3_en}

    .. image:: {volcano_4_en}

    .. image:: {volcano_5_en}

        """.format( unique_reads_en=unique_reads_en, rRNA_metrics_en=rRNA_metrics_en, read_distrib_en=read_distrib_en
        , genebody_hm_en=genebody_hm_en, genebody_cv_en=genebody_cv_en, heatmapSF_en=heatmapSF_en, heatmapSS_en=heatmapSS_en,
        pca_1_en=pca_1_en,pca_2_en=pca_2_en,pca_3_en=pca_3_en,heatmapSSC_en=heatmapSSC_en,volcano_1_en=volcano_1_en,volcano_2_en=volcano_2_en,volcano_3_en=volcano_3_en,volcano_4_en=volcano_4_en,volcano_5_en=volcano_5_en ), output[0], metadata="Molecular Biology Core Facilities, DFCI", **input)


rule run_STAR:
    input:
        get_fastq
    output:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        counts="analysis/STAR/{sample}/{sample}.counts.tab",
        log_file="analysis/STAR/{sample}/{sample}.Log.final.out"
    params:
        stranded=strand_command,
        prefix=lambda wildcards: "analysis/STAR/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_index]}"
	" --sjdbGTFfile {config[gtf_file]}"
        " --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.prefix}."
	"  --outSAMstrandField intronMotif"
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000 --quantMode GeneCounts"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"
        " && /usr/bin/samtools index {output.bam}"


rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR/{sample}/{sample}.Log.final.out", sample=ordered_sample_list ),
        star_gene_count_files=expand( "analysis/STAR/{sample}/{sample}.counts.tab", sample=ordered_sample_list )
    output:
        csv="analysis/STAR/STAR_Align_Report.csv",
        png="analysis/STAR/STAR_Align_Report.png",
        gene_counts="analysis/STAR/STAR_Gene_Counts.csv"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl snakemake/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript snakemake/scripts/map_stats.R {output.csv} {output.png}" )
        shell( "perl snakemake/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

rule run_cufflinks:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking"
    threads: 4
    params:
        library_command=cuff_command
    shell:
        "cufflinks -o analysis/cufflinks/{wildcards.sample} -p {threads} -G {config[gtf_file]} {params.library_command} {input}"
        " && mv analysis/cufflinks/{wildcards.sample}/genes.fpkm_tracking {output}"
        " && mv analysis/cufflinks/{wildcards.sample}/isoforms.fpkm_tracking analysis/cufflinks/{wildcards.sample}/{wildcards.sample}.isoforms.fpkm_tracking"

rule generate_cuff_matrix:
    input:
        cuff_gene_fpkms=expand( "analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking", sample=ordered_sample_list )
    output:
        "analysis/cufflinks/Cuff_Gene_Counts.csv"
    run:
        fpkm_files= " -f ".join( input.cuff_gene_fpkms )
        shell( "perl snakemake/scripts/raw_and_fpkm_count_matrix.pl -c -f {fpkm_files} 1>{output}" )


rule run_STAR_fusion:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam" #just to make sure STAR output is available before STAR_Fusion
    output:
        "analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final"
    log:
        "analysis/STAR_Fusion/{sample}/{sample}.star_fusion.log"
    shell:
        "STAR-Fusion --chimeric_junction analysis/STAR/{wildcards.sample}/{wildcards.sample}.Chimeric.out.junction "
        "--genome_lib_dir {config[genome_lib_dir]} --output_dir analysis/STAR_Fusion/{wildcards.sample} >& {log}"
        " && mv analysis/STAR_Fusion/{wildcards.sample}/star-fusion.fusion_candidates.final {output}"
        " && mv analysis/STAR_Fusion/{wildcards.sample}/star-fusion.fusion_candidates.final.abridged"
        " analysis/STAR_Fusion/{wildcards.sample}/{wildcards.sample}.fusion_candidates.final.abridged"


rule read_distrib_qc:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/read_distrib/{sample}.txt"
    shell:
        "{config[python2]} {config[rseqc_path]}/read_distribution.py"
        " --input-file={input}"
        " --refgene={config[bed_file]} 1>{output}"

rule read_distrib_qc_matrix:
    input:
        read_distrib_files=expand( "analysis/RSeQC/read_distrib/{sample}.txt", sample=ordered_sample_list )
    output:
        matrix="analysis/RSeQC/read_distrib/read_distrib.matrix.tab",
        png="analysis/RSeQC/read_distrib/read_distrib.png"
    run:
        file_list_with_flag = " -f ".join( input.read_distrib_files )
        shell( "perl snakemake/scripts/read_distrib_matrix.pl -f {file_list_with_flag} 1>{output.matrix}" )
        shell( "Rscript snakemake/scripts/read_distrib.R {output.matrix} {output.png}" )


rule down_sample:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/gene_body_cvg/downsample/{sample}.downsample.sorted.bam"
    shell:
        "java -jar {config[picard_path]}/DownsampleSam.jar INPUT={input} OUTPUT={output}"
        " PROBABILITY=0.1"
        " && samtools index {input}"

rule gene_body_cvg_qc:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png",
        "analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r"
    shell:
        "{config[python2]} {config[rseqc_path]}/geneBody_coverage.py -i {input} -r {config[bed_file]}"
        " -f png -o analysis/RSeQC/gene_body_cvg/{wildcards.sample}/{wildcards.sample}"


rule plot_gene_body_cvg:
    input:
        expand("analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", sample=ordered_sample_list )
    output:
        rscript="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.r",
        png="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        png_curves="analysis/RSeQC/gene_body_cvg/geneBodyCoverage.curves.png"
    shell:
        "perl snakemake/scripts/plot_gene_body_cvg.pl --rfile {output.rscript} --png {output.png} --curves_png {output.png_curves} {input}"
        " && Rscript {output.rscript}"

rule junction_saturation:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf"
    shell:
        "{config[python2]} {config[rseqc_path]}/junction_saturation.py -i {input} -r {config[bed_file]}"
        " -o analysis/RSeQC/junction_saturation/{wildcards.sample}/{wildcards.sample}"


rule collect_insert_size:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/insert_size/{sample}/{sample}.histogram.pdf"
    shell:
        "java -jar {config[picard_path]}/CollectInsertSizeMetrics.jar"
        " H={output} I={input} O=analysis/RSeQC/insert_size/{wildcards.sample}/{wildcards.sample} R={config[ref_fasta]}"


rule run_rRNA_STAR:
    input:
        get_fastq
    output:
        bam="analysis/STAR_rRNA/{sample}/{sample}.sorted.bam",
        log_file="analysis/STAR_rRNA/{sample}/{sample}.Log.final.out"
    params:
        stranded=rRNA_strand_command,
        prefix=lambda wildcards: "analysis/STAR_rRNA/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    shell:
        "STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_rRNA_index]}"
        " --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.prefix}."
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && samtools index {output.bam}"


rule generate_rRNA_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR_rRNA/{sample}/{sample}.Log.final.out", sample=ordered_sample_list )
    output:
        csv="analysis/STAR_rRNA/STAR_rRNA_Align_Report.csv",
        png="analysis/STAR_rRNA/STAR_rRNA_Align_Report.png"
    run:
        log_files = " -l ".join( input.star_log_files )
        shell( "perl snakemake/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript snakemake/scripts/map_stats_rRNA.R {output.csv} {output.png}" )

rule get_chrom_size:
    output:
        "analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
    params:
        config['reference']
    shell:
        "fetchChromSizes {params} 1>{output}"

#MAHESH's
rule bam_to_bigwig:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        chrom_size="analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
    output:
        "analysis/bam2bw/{sample}/{sample}.bw"
    params:
        "analysis/bam2bw/{sample}/{sample}"
    shell:
        "bedtools genomecov -bg -split -ibam {input.bam} -g {input.chrom_size} 1> {params}.bg"
        " && bedSort {params}.bg {params}.sorted.bg"
        " && bedGraphToBigWig {params}.sorted.bg {input.chrom_size} {output}"

#LEN:
# rule bam_to_bg:
#     input:
#         bam="analysis/STAR/{sample}/{sample}.sorted.bam",
#         chrom_size="analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
#     output:
#         "analysis/bam2bw/{sample}/{sample}.bg"
#     shell:
#         "bedtools genomecov -bg -split -ibam {input.bam} -g {input.chrom_size} 1> {output}"

# rule sort_bg:
#     input:
#         "analysis/bam2bw/{sample}/{sample}.bg"
#     output:
#         "analysis/bam2bw/{sample}/{sample}.sorted.bg"
#     shell:
#         "bedSort {input} {output}"

# rule bg_to_bw:
#     input:
#         bg="analysis/bam2bw/{sample}/{sample}.sorted.bg",
#         chrom_size="analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
#     output:
#         "analysis/bam2bw/{sample}/{sample}.bw"
#     shell:
#         "bedGraphToBigWig {input.bg} {input.chrom_size} {output}"

#LEN:
rule pca_plot:
    input:
        rpkmFile="analysis/cufflinks/Cuff_Gene_Counts.csv",
        annotFile=config['metasheet']
    output:
        pca_plot_out="analysis/plots/pca_plot.pdf",
        png_dir="analysis/plots/images/"
#    shell:
#        "scripts/pca_plot.R"
    run:
        shell("Rscript snakemake/scripts/pca_plot.R {input.rpkmFile} {input.annotFile} {output.pca_plot_out} {output.png_dir}")

rule heatmapSS_plot:
    input:
        rpkmFile="analysis/cufflinks/Cuff_Gene_Counts.csv",
        annotFile=config['metasheet']
    output:
        ss_plot_out="analysis/plots/heatmapSS_plot.pdf",
        ss_txt_out="analysis/plots/heatmapSS.txt"
    run:
        shell("Rscript snakemake/scripts/heatmapSS_plot.R {input.rpkmFile} {input.annotFile} {output.ss_plot_out} {output.ss_txt_out}")

rule heatmapSF_plot:
    input:
        rpkmFile="analysis/cufflinks/Cuff_Gene_Counts.csv",
        annotFile=config['metasheet']
    output:
        sf_plot_out="analysis/plots/heatmapSF_plot.pdf",
        sf_txt_out="analysis/plots/heatmapSF.txt"
    run:
        shell("Rscript snakemake/scripts/heatmapSF_plot.R {input.rpkmFile} {input.annotFile} {output.sf_plot_out} {output.sf_txt_out}")

#PART 2.2- diffexp w/ DEseq
#based on tosh's coppRhead/Snakefile

## Get the raw counts for each sample
rule gfold_count:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        reference_genes=config["gtf_file"],
    output:
        "analysis/gfold/{sample}/{sample}.read_cnt.txt"
    shell:
        "samtools view {input.bam} | "
        "gfold count -ann {input.reference_genes} -tag stdin -o {output}"

## Extract comparisons from the metadata file and perform gfold diff
def get_column(comparison):
    return metadata["comp_{}".format(comparison)]

def get_comparison(name, group):
    comp = get_column(name)
    return metadata[comp == group].index

def get_samples(wildcards):
    comp = get_column(wildcards.comparison)
    return comp.dropna().index

## Perform Limma and DEseq on comparisons
rule limma_and_deseq:
    input:
        counts = lambda wildcards: expand("analysis/gfold/{sample}/{sample}.read_cnt.txt", sample=get_samples(wildcards))
    output:
        limma = "analysis/diffexp/{comparison}/{comparison}.limma.txt",
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.txt",
        deseqSum = "analysis/diffexp/{comparison}/{comparison}.deseq.sum.csv"
    params:
        s1=lambda wildcards: ",".join(get_comparison(wildcards.comparison, 1)),
        s2=lambda wildcards: ",".join(get_comparison(wildcards.comparison, 2))
#    script:
#        "scripts/DEseq.R"

    run:
        counts = " ".join(input.counts)
        shell("Rscript snakemake/scripts/DEseq.R \"{counts}\" \"{params.s1}\" \"{params.s2}\" {output.limma} {output.deseq} {output.deseqSum}")

rule fetch_DE_gene_list:
    input:
        deseq_file_list=expand("analysis/diffexp/{comparison}/{comparison}.deseq.txt",comparison=comparisons)
    output:
        "analysis/diffexp/de_summary.csv"
    run:
        deseq_file_string = ' -f '.join(input.deseq_file_list)
        shell("perl snakemake/scripts/get_de_summary_table.pl -f {deseq_file_string} 1>{output}") 

#Generate volcano plots for each comparison
rule volcano_plot:
    input:
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.txt",
    output:
        plot = "analysis/diffexp/{comparison}/{comparison}_volcano.pdf",
        png = "analysis/plots/images/{comparison}_volcano.png"
    run:
        shell("Rscript snakemake/scripts/volcano_plot.R {input.deseq} {output.plot} {output.png}")
