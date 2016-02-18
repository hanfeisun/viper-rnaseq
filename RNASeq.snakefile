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
        expand( "analysis/bam2bw/{sample}/{sample}.bw", sample=ordered_sample_list ),
        expand( "analysis/gfold/{sample}/{sample}.read_cnt.txt", sample=ordered_sample_list ),
        expand("analysis/diffexp/{comparison}/{comparison}.deseq.txt", comparison=comparisons),
        expand("analysis/diffexp/{comparison}/{comparison}_volcano.pdf", comparison=comparisons),
        "analysis/diffexp/de_summary.png",
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
        heatmapSF_plot="analysis/plots/images/heatmapSF_plot.png",
        heatmapSS_plot="analysis/plots/images/heatmapSS_plot.png",
       # heatmapSS_cluster="analysis/plots/images/heatmapSS_cluster.png",
        DEsummary_plot="analysis/diffexp/de_summary.png",
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
        heatmapSF_en = data_uri( input.heatmapSF_plot )
        heatmapSS_en = data_uri( input.heatmapSS_plot )
        #heatmapSSC_en = data_uri(input.heatmapSS_cluster)
        DEsummary_en = data_uri (input.DEsummary_plot)

        ### Getting all pdf reports ###
        import os
        import re
        import glob
        pdf_list = {}
        ignore_pdf = re.compile('.*RSeQC/junction_saturation.*')
        for root,dirs,files in os.walk('./analysis'):
            for file in files:
                if( (not bool(ignore_pdf.match(os.path.join(root,file)))) and (not file.startswith('.')) and (file.endswith('.pdf')) ):
                    pdf_list[file[:-4]] = os.path.join(root,file)

        ### Getting all pca and volcano plots
        pca_png_list = []
        volcano_list = []

        for pca_plot in sorted(glob.glob("./analysis/plots/images/pca_plot*.png")):
            pca_png_list.append(data_uri(pca_plot))

        for volcano_plot in glob.glob("./analysis/plots/images/*_volcano.png"):
            volcano_list.append(data_uri(volcano_plot))

        pca_string = "\n\t.. image:: " + "\n\t.. image:: ".join(pca_png_list[:-1]) + "\n"
        pca_var_string = "\n\t.. image:: " + pca_png_list[-1] + "\n"
        volcano_string = "\n\t.. image:: " + "\n\t.. image:: ".join(volcano_list) + "\n"



        report("""
==============================================
VIPER: Visualization Pipeline for RNAseq
==============================================


Alignment Summary
=================
    Raw reads were mapped or aligned to the reference organism using `STAR software`_.

    .. _STAR software: https://github.com/alexdobin/STAR


    The **uniquely mapped read counts** and the **total read counts** for all the samples are summarized in the following image. In most cases, more than 70% of the reads should uniquely map to the reference genome.
    Contamination, low quality sequencing data, or poor library contruction may result in <60% uniquely aligned reads.

    .. image:: {unique_reads_en}

Library Prep Quality Metrics
=============================
Read Distribution QC
^^^^^^^^^^^^^^^^^^^^^^^^^^
    This graph displays the disibution of reads mapped to **features** across the genome for each sample. Distribution profiles should be similar across samples.
    A sample with a siginficantly different distribution profile may indicate that it was initially a poor sample or that there was a problem during library preparation.
    **mRNAseq libraries will typically contain less than 20% intronic mapped reads** whereas **totalRNAseq libraries may have greater than 40% intronic mapped reads**.

    .. image:: {read_distrib_en}

rRNA removal QC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    This graph displays the percentage of reads mapping to ribosomal RNA reference sequences.
    Most RNAseq library prep methods are designed to avoid sampling ribosomal RNAs that in some cases can greater than 80% of total RNA.
    *mRNAseq* library prep methods address this issue by using oligo-dT RT priming or capture beads to enrich for polyadenylated transcripts. *totalRNAseq* library prep methods address this problem by selectivel depleting rRNAs.
    Analyzing reads mapping to rRNA is a useful metric to assess the effectiveness of either mRNA enrichment or rRNA depletion.
    For **mRNAseq prep methods ~5% rRNA** read mapping is typical. It should be **<10% for totalRNAseq**. Higher percentages may indicate poor mRNA enrichment or rRNA depletion.

    .. image:: {rRNA_metrics_en}


Genebody Coverage
^^^^^^^^^^^^^^^^^
    For accurate gene expression quantification, mapped reads should be evenly distributed across genebodies.
    Significantly skewed profiles (5' or 3') may introduce quantification bias and/or represent poor quality library preparation.\n
    For example, mRNAseq library preps typically use oligo-dT beads to capture mature transcripts and can be prone to 3' bias in genebody coverage if degraded RNA \(RIN < 7\) is used as input. This may result in inaccurate gene quantification and the following graphs will help diagnose.
    There are other prep methods that may result in 5' bias too. Ideally, coverage should be uniform across the genebody. The line plots should look like this: "∩"

    Figures generated using `RSeQC software`_.

    .. _RSeQC software: http://rseqc.sourceforge.net

    **Line Plot**

    .. image:: {genebody_cv_en}

    **Heatmap**\n
    This graphic may facilitate identification of biased samples.\n
    Scale: Blue = 0 Pink =1

    .. image:: {genebody_hm_en}


Experimental Quality Control
===============================

Principle Component Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    High dimensional expression data are mathmatically reduced to principle
    components that can be used to describe variation across samples in fewer dimensions to allow human interpretation.
    Principle component 1 \(PC1\) accounts for the most amount of variation across samples, PC2 the second most, and so on. These PC1 vs PC2 plots
    are colored by sample annotation to demontrate how samples cluster together \(or not\) in reduced dimensional space.
    For more detailed description of Princilple Component Analysis, start with `wikipedia`_.

    .. _wikipedia: https://en.wikipedia.org/wiki/Principal_component_analysis

{pca_string}

    This plot indicates how much of the overall variance is described by the principle components in descending order.

{pca_var_string}


Sample-to-Sample Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    This heatmap displays hierarchical clustering of spearman rank correlations across samples.
    Highly correlated samples will cluster together and be represented in red. Samples that do not correlate will be represented in blue.

    .. image:: {heatmapSS_en}

    Sample-to-Sample data matrix is /analysis/plots/heatmapSS.txt

Sample-Feature Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    This heatmap illustrates hierarchical clustering of samples based on the top 5 percent or roughly 1000 most variable genes or "features."
    The top colomn color annotations are presented to help identify how sample types, groups, or treatments are clustering together \(or not\).


    .. image:: {heatmapSF_en}

    What are *these* genes?
    Data used to generate this sample-feature graphic are in /analysis/plots/heatmapSF.txt

Differential Gene expression
============================
    Differential gene expression analysis was performed using both `limma`_ and `DESeq2`_.\n
    Full analysis output tables are are available in /analysis/diffexp/comparison_of_interest

    .. _limma: https://www.bioconductor.org/packages/3.3/bioc/vignettes/limma/inst/doc/usersguide.pdf

    .. _DESeq2: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/

    This summary image shows the number of genes that are up regulated and down regulated across each comparison at different adjusted P-value cut-offs.

    .. image:: {DEsummary_en}

Volcano Plots
^^^^^^^^^^^^^^
    Volcano plots are commonly used graphical representations of differentially expressed genes and statistical significance.
    These scatter plots show log2 fold change versus P-value for all genes detected. Each data point represents a gene. Genes are colored red if the log2 fold change is greater than one \(log2fc > 1\). Genes are colored blue if the log2 fold change is less than negative one \(log2fc < -1\).
    The plot title indicates the direction of the comparison.
    For example, "treatment_vs_control" indicates that genes colored red are up-regulated in the treatment condition compared to the control condition with a statistically significant P-value.

{volcano_string}

        """.format( unique_reads_en=unique_reads_en, rRNA_metrics_en=rRNA_metrics_en, read_distrib_en=read_distrib_en
        , genebody_hm_en=genebody_hm_en, genebody_cv_en=genebody_cv_en, heatmapSF_en=heatmapSF_en, heatmapSS_en=heatmapSS_en,
        DEsummary_en=DEsummary_en,pca_string=pca_string,pca_var_string=pca_var_string,volcano_string=volcano_string ),
        output[0], metadata="Molecular Biology Core Facilities, DFCI", **{'Copyrights:':"./snakemake/mbcf.jpg"})


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
        " && touch {output}" # For some sample, final.abridged is created but not .final file; temp hack before further investigate into this


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
        " && if [ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/ERCC/input/ERCC92.chromInfo ]; then cat /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/ERCC/input/ERCC92.chromInfo 1>>{output}; fi"

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
    params:
        RPKM_threshold = config["RPKM_threshold"],
        min_num_samples_expressing_at_threshold = config["min_num_samples_expressing_at_threshold"],
        SSnumgenes = config["SSnumgenes"]
    run:
        shell("Rscript snakemake/scripts/heatmapSS_plot.R {input.rpkmFile} {input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} {params.SSnumgenes} {output.ss_plot_out} {output.ss_txt_out}")

rule heatmapSF_plot:
    input:
        rpkmFile="analysis/cufflinks/Cuff_Gene_Counts.csv",
        annotFile=config['metasheet']
    output:
        sf_plot_out="analysis/plots/heatmapSF_plot.pdf",
        sf_txt_out="analysis/plots/heatmapSF.txt"
    params:
        RPKM_threshold = config["RPKM_threshold"],
        min_num_samples_expressing_at_threshold = config["min_num_samples_expressing_at_threshold"],
        SFnumgenes = config["SFnumgenes"],
        num_kmeans_clust = config["num_kmeans_clust"]
    run:
        shell("Rscript snakemake/scripts/heatmapSF_plot.R {input.rpkmFile} {input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} {params.SFnumgenes} {params.num_kmeans_clust} {output.sf_plot_out} {output.sf_txt_out}")

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
        csv="analysis/diffexp/de_summary.csv",
        png="analysis/diffexp/de_summary.png"
    run:
        deseq_file_string = ' -f '.join(input.deseq_file_list)
        shell("perl snakemake/scripts/get_de_summary_table.pl -f {deseq_file_string} 1>{output.csv}")
        shell("Rscript snakemake/scripts/de_summary.R {output.csv} {output.png}")

#Generate volcano plots for each comparison
rule volcano_plot:
    input:
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.txt",
    output:
        plot = "analysis/diffexp/{comparison}/{comparison}_volcano.pdf",
        png = "analysis/plots/images/{comparison}_volcano.png"
    run:
        shell("Rscript snakemake/scripts/volcano_plot.R {input.deseq} {output.plot} {output.png}")
