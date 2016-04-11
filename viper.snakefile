# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
from collections import defaultdict
import pandas as pd
import yaml

#-----     CONFIG SET UP    ------#
configfile: "config.yaml"

with open("ref.yaml","r") as ref_file:
    ref_info = yaml.safe_load(ref_file) 

for k,v in ref_info.items():
    config[k] = v

config["samples"] = config["the_samples"]

for k in ["RPKM_threshold","min_num_samples_expressing_at_threshold","SSnumgenes","SFnumgenes","num_kmeans_clust","filter_mirna","snp_scan_genome"]:
    config[k] = str(config[k])
#----   END OF CONFIG SET UP -----#

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


metadata = pd.read_table(config['metasheet'], index_col=0, sep=',')
comparisons = comparison=[c[5:] for c in metadata.columns if c.startswith("comp_")]

metacols = [c for c in metadata.columns if c.lower()[:4] != 'comp']

# Mahesh changing the metasheet match with config info to using pandas #
file_info = { sampleName : config["samples"][sampleName] for sampleName in metadata.index }
ordered_sample_list = metadata.index
run_fusion= True if len(config["samples"][metadata.index[0]]) == 2 else False

if( run_fusion ):
    if( config["stranded"] ):
        strand_command = " --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000"
    else:
        strand_command = " --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif"

#GENERATE snp regions list:
snp_regions = ['chr6', 'genome'] if ('snp_scan_genome' in config) and config['snp_scan_genome'] == 'true' else ['chr6']

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


rule target:
    input:
        expand( "analysis/cufflinks/{K}/{K}.genes.fpkm_tracking", K=ordered_sample_list ),
        "analysis/STAR/STAR_Align_Report.csv",
        "analysis/STAR/STAR_Align_Report.png",
        "analysis/STAR/STAR_Gene_Counts.csv",
        "analysis/cufflinks/Cuff_Gene_Counts.csv",
        "analysis/plots/pca_plot.pdf",
        expand("analysis/plots/images/pca_plot_{metacol}.png", metacol=metacols),
        "analysis/plots/heatmapSS_plot.pdf",
        "analysis/plots/heatmapSF_plot.pdf",
        expand( "analysis/RSeQC/read_distrib/{sample}.txt", sample=ordered_sample_list ),
        "analysis/RSeQC/read_distrib/read_distrib.png",
        expand( "analysis/RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png", sample=ordered_sample_list ),
        "analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        expand( "analysis/RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf", sample=ordered_sample_list ),
        expand( "analysis/bam2bw/{sample}/{sample}.bw", sample=ordered_sample_list ),
        expand("analysis/diffexp/{comparison}/{comparison}.deseq.csv", comparison=comparisons),
        expand("analysis/diffexp/{comparison}/{comparison}_volcano.pdf", comparison=comparisons),
        "analysis/diffexp/de_summary.png",
        expand( "analysis/snp/{sample}/{sample}.snp.{region}.txt", sample=ordered_sample_list, region=snp_regions ),
        expand( "analysis/snp/snp_corr.{region}.txt", region=snp_regions ),
        expand( "analysis/plots/sampleSNPcorr_plot.{region}.png", region=snp_regions),
        fusion_output,
        insert_size_output,
        rRNA_metrics
    message: "Compiling all output"
        
#["analysis/plots/correlation_plot.pdf", "analysis/plots/correlation_table.csv", "analysis/plots/upvenn_plot.pdf", "analysis/plots/downvenn_plot.pdf"] if len(comparisons) >= 2 else []


def get_sphinx_report():
    import os
    import glob
    from snakemake.report import data_uri
    file_dict = {
        'align_report': "analysis/STAR/STAR_Align_Report.png",
        'rRNA_report': "analysis/STAR_rRNA/STAR_rRNA_Align_Report.png",
        'read_distrib': "analysis/RSeQC/read_distrib/read_distrib.png",
        'gb_cov_heatmap': "analysis/RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        'gb_cov_curves': "analysis/RSeQC/gene_body_cvg/geneBodyCoverage.curves.png",
        'heatmapSF_plot': "analysis/plots/images/heatmapSF_plot.png",
        'heatmapSS_plot': "analysis/plots/images/heatmapSS_plot.png",
        'heatmapSS_cluster': "analysis/plots/images/heatmapSS_cluster.png",
        'DEsummary_plot': "analysis/diffexp/de_summary.png",
    }
    copy_file_dict = {}
    for key in file_dict.keys():
        copy_file_dict[key] = file_dict[key]
    for file_token in file_dict.keys():
        if not os.path.isfile(file_dict[file_token]):
            del copy_file_dict[file_token]
        else:
            copy_file_dict[file_token] = data_uri(copy_file_dict[file_token])
    file_dict = copy_file_dict
    pca_png_list = []
    volcano_list = []
    SF_png_list = []

    for pca_plot in sorted(glob.glob("./analysis/plots/images/pca_plot*.png")):
        pca_png_list.append(data_uri(pca_plot))

    for volcano_plot in glob.glob("./analysis/plots/images/*_volcano.png"):
        volcano_list.append(data_uri(volcano_plot))

    for SF_plot in sorted(glob.glob("./analysis/plots/images/heatmapSF_*_plot.png")):
        SF_png_list.append(data_uri(SF_plot))    

    if pca_png_list:
        file_dict['pca_png_list'] = pca_png_list
    if volcano_list:
        file_dict['volcano_png_list'] = volcano_list
    if SF_png_list:
        file_dict['sf_png_list'] = SF_png_list

    report = """
==============================================
VIPER: Visualization Pipeline for RNAseq
==============================================


Alignment Summary
=================
    Raw reads were mapped or aligned to the reference organism using `STAR software`_.

    .. _STAR software: https://github.com/alexdobin/STAR


    The **uniquely mapped read counts** and the **total read counts** for all the samples are summarized in the following image. In most cases, more than 70% of the reads should uniquely map to the reference genome.
    Contamination, low quality sequencing data, or poor library contruction may result in <60% uniquely aligned reads.

"""
    if 'align_report' in file_dict:
        report += "\n\t.. image:: " + file_dict['align_report'] + "\n";

    report += "\n"
    report += """
Library Prep Quality Metrics
=============================
Read Distribution QC
^^^^^^^^^^^^^^^^^^^^^^^^^^
    This graph displays the disibution of reads mapped to **features** across the genome for each sample. Distribution profiles should be similar across samples.
    A sample with a siginficantly different distribution profile may indicate that it was initially a poor sample or that there was a problem during library preparation.
    **mRNAseq libraries will typically contain less than 20% intronic mapped reads** whereas **totalRNAseq libraries may have greater than 40% intronic mapped reads**.
"""

    if 'read_distrib' in file_dict:
        report += "\n\t.. image:: " + file_dict['read_distrib'] + "\n";

    report += "\n"
    report += """
rRNA removal QC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    This graph displays the percentage of reads mapping to ribosomal RNA reference sequences. Most RNAseq library prep methods are designed to avoid sampling ribosomal RNAs which typically represent greater than 80% of total RNA. If rRNA removal was effective, less than 5% of the reads should map to rRNA sequences and for mRNA libraries fewer than 1%.
"""

    if 'rRNA_report' in file_dict:
        report += "\n\t.. image:: " + file_dict['rRNA_report'] + "\n";

    report += "\n"
    report += """
Genebody Coverage
^^^^^^^^^^^^^^^^^
    For accurate gene expression quantification, mapped reads should be evenly distributed across genebodies.
    Significantly skewed profiles (5' or 3') may introduce quantification bias and/or represent poor quality library preparation.\n
    For example, mRNAseq library preps typically use oligo-dT beads to capture mature transcripts and can be prone to 3' bias in genebody coverage if degraded RNA \(RIN < 7\) is used as input. This may result in inaccurate gene quantification and the following graphs will help diagnose.
    There are other prep methods that may result in 5' bias too. Ideally, coverage should be uniform across the genebody. The line plots should look like this: "∩"

    Figures generated using `RSeQC software`_.

    .. _RSeQC software: http://rseqc.sourceforge.net

    **Line Plot**
"""

    if 'gb_cov_curves' in file_dict:
        report += "\n\t.. image:: " + file_dict['gb_cov_curves'] + "\n";

    report += "\n"
    report += """
    **Heatmap**\n
    This graphic may facilitate identification of biased samples.\n
    Scale: Blue = 0 Pink =1

"""

    if 'gb_cov_heatmap' in file_dict:
        report += "\n\t.. image:: " + file_dict['gb_cov_heatmap'] + "\n";

    report += "\n"
    report += """
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

"""

    if 'pca_png_list' in file_dict:
        if len(file_dict['pca_png_list']) > 1:
            report += "\n\t.. image:: " + "\n\t.. image:: ".join(file_dict['pca_png_list'][:-1]) + "\n"
            report += "\n\t" + 'This plot indicates how much of the overall variance is described by the principle components in descending order.' + "\n\n\t.. image:: " + file_dict['pca_png_list'][-1] + "\n"
        else:
            report += "\n\t.. image:: " + "\n\t.. image:: ".join(file_dict['pca_png_list'][0]) + "\n"

    report += "\n"
    report += """
Sample-to-Sample Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    This heatmap displays hierarchical clustering of spearman rank correlations across samples.
    Highly correlated samples will cluster together and be represented in red. Samples that do not correlate will be represented in blue.

"""

    if 'heatmapSS_plot' in file_dict:
        report += "\n\n\t.. image:: " + file_dict['heatmapSS_plot'] + "\n\n\t" + 'Sample-to-Sample data matrix is /analysis/plots/heatmapSS.txt' + "\n"

    report += "\n"
    report += """
Sample-Feature Correlation Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    This heatmap illustrates hierarchical clustering of samples based on the top 5 percent or roughly 1000 most variable genes or "features."
    The top colomn color annotations are presented to help identify how sample types, groups, or treatments are clustering together \(or not\).


"""

    if 'heatmapSF_plot' in file_dict:
        report += "\n\n\t.. image:: " + file_dict['heatmapSF_plot'] + "\n";

    if 'sf_png_list' in file_dict:
        report += "\n\t.. image:: " + "\n\t.. image:: ".join(file_dict['sf_png_list'][:]) + "\n"

    report += "\n\n\t" + 'What are *these* genes?' + "\n\n\t" + 'Data used to generate this sample-feature graphic are in /analysis/plots/heatmapSF.txt' + "\n"

    report += "\n"

    report += """                
Differential Gene expression
============================
    Differential gene expression analysis was performed using both `limma`_ and `DESeq2`_.\n
    Full analysis output tables are are available in /analysis/diffexp/comparison_of_interest

    .. _limma: https://www.bioconductor.org/packages/3.3/bioc/vignettes/limma/inst/doc/usersguide.pdf

    .. _DESeq2: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049/

    This summary image shows the number of genes that are up regulated and down regulated across each comparison at different adjusted P-value cut-offs.

"""

    if 'DEsummary_plot' in file_dict:
        report += "\n\n\t.. image:: " + file_dict['DEsummary_plot'] + "\n"

    report += "\n"
    report += """
Volcano Plots
^^^^^^^^^^^^^^
    Volcano plots are commonly used graphical representations of differentially expressed genes and statistical significance.
    These scatter plots show log2 fold change versus P-value for all genes detected. Each data point represents a gene. Genes are colored red if the log2 fold change is greater than one \(log2fc > 1\). Genes are colored blue if the log2 fold change is less than negative one \(log2fc < -1\).
    The plot title indicates the direction of the comparison.
    For example, "treatment_vs_control" indicates that genes colored red are up-regulated in the treatment condition compared to the control condition with a statistically significant P-value.

"""

    if 'volcano_png_list' in file_dict:
        report += "\n\n\t.. image:: " + "\n\n\t.. image:: ".join(file_dict['volcano_png_list'][:]) + "\n"

    return report + "\n"

rule generate_report:
    input:
    output:
        "report.html"
    message: "Generating VIPER report"
    run:
        from snakemake.utils import report
        sphinx_str = get_sphinx_report()
        report(sphinx_str, output[0], metadata="Molecular Biology Core Facilities, DFCI", **{'Copyrights:':"./viper/mbcf.jpg"})


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
    message: "Running STAR Alignment on {wildcards.sample}"
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
    message: "Generating STAR report"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl viper/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript viper/scripts/map_stats.R {output.csv} {output.png}" )
        shell( "perl viper/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

rule run_cufflinks:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/cufflinks/{sample}/{sample}.genes.fpkm_tracking"
    threads: 4
    message: "Running Cufflinks on {wildcards.sample}"
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
    message: "Generating expression matrix using cufflinks counts"
    run:
        fpkm_files= " -f ".join( input.cuff_gene_fpkms )
        shell( "perl viper/scripts/raw_and_fpkm_count_matrix.pl -c -f {fpkm_files} 1>{output}" )


rule run_STAR_fusion:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam" #just to make sure STAR output is available before STAR_Fusion
    output:
        "analysis/STAR_Fusion/{sample}/{sample}.fusion_candidates.final"
    log:
        "analysis/STAR_Fusion/{sample}/{sample}.star_fusion.log"
    message: "Running STAR fusion on {wildcards.sample}"
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
    message: "Running RseQC read distribution on {wildcards.sample}"
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
    message: "Creating RseQC read distribution matrix"
    run:
        file_list_with_flag = " -f ".join( input.read_distrib_files )
        shell( "perl viper/scripts/read_distrib_matrix.pl -f {file_list_with_flag} 1>{output.matrix}" )
        shell( "Rscript viper/scripts/read_distrib.R {output.matrix} {output.png}" )

rule down_sample:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/gene_body_cvg/downsample/{sample}.downsample.sorted.bam"
    message: "Running RseQC downsample gene body coverage for {wildcards.sample}"
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
    message: "Creating gene body coverage curves"
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
    message: "Plotting gene body coverage"
    shell:
        "perl viper/scripts/plot_gene_body_cvg.pl --rfile {output.rscript} --png {output.png} --curves_png {output.png_curves} {input}"
        " && Rscript {output.rscript}"

rule junction_saturation:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf"
    message: "Determining junction saturation for {wildcards.sample}"
    shell:
        "{config[python2]} {config[rseqc_path]}/junction_saturation.py -i {input} -r {config[bed_file]}"
        " -o analysis/RSeQC/junction_saturation/{wildcards.sample}/{wildcards.sample}"


rule collect_insert_size:
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/RSeQC/insert_size/{sample}/{sample}.histogram.pdf"
    message: "Collecting insert size for {wildcards.sample}"
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
    message: "Running rRNA STAR for {wildcards.sample}"
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
    message: "Generating STAR rRNA report"
    run:
        log_files = " -l ".join( input.star_log_files )
        shell( "perl viper/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript viper/scripts/map_stats_rRNA.R {output.csv} {output.png}" )

rule get_chrom_size:
    output:
        "analysis/bam2bw/" + config['reference'] + ".Chromsizes.txt"
    params:
        config['reference']
    message: "Fetching chromosome sizes"
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
    message: "Converting {wildcards.sample} bam to bigwig"
    shell:
        "bedtools genomecov -bg -split -ibam {input.bam} -g {input.chrom_size} 1> {params}.bg"
        " && bedSort {params}.bg {params}.sorted.bg"
        " && bedGraphToBigWig {params}.sorted.bg {input.chrom_size} {output}"

rule pca_plot:
    input:
        rpkmFile="analysis/cufflinks/Cuff_Gene_Counts.csv",
        annotFile=config['metasheet']
    output:
        expand("analysis/plots/images/pca_plot_{metacol}.png", metacol=metacols),
        pca_plot_out="analysis/plots/pca_plot.pdf"
    params:
        RPKM_threshold = config["RPKM_threshold"],
        min_num_samples_expressing_at_threshold = config["min_num_samples_expressing_at_threshold"],
        filter_mirna = config["filter_mirna"],
        SSnumgenes = config["SSnumgenes"]
    message: "Generating PCA plots"
#    shell:
#        "scripts/pca_plot.R"
    run:
        shell("Rscript viper/scripts/pca_plot.R {input.rpkmFile} {input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} {params.filter_mirna} {params.SSnumgenes} {output.pca_plot_out}")

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
        filter_mirna = config["filter_mirna"],
        SSnumgenes = config["SSnumgenes"]
    message: "Generating Sample-Sample Heatmap"
    run:
        shell("Rscript viper/scripts/heatmapSS_plot.R {input.rpkmFile} {input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} {params.filter_mirna} {params.SSnumgenes} {output.ss_plot_out} {output.ss_txt_out}")

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
        filter_mirna = config["filter_mirna"],
        SFnumgenes = config["SFnumgenes"],
        num_kmeans_clust = config["num_kmeans_clust"]
    message: "Generating Sample-Feature heatmap"
    run:
        shell("Rscript viper/scripts/heatmapSF_plot.R {input.rpkmFile} {input.annotFile} {params.RPKM_threshold} {params.min_num_samples_expressing_at_threshold} {params.filter_mirna} {params.SFnumgenes} {params.num_kmeans_clust} {output.sf_plot_out} {output.sf_txt_out}")

#PART 2.2- diffexp w/ DEseq
#based on tosh's coppRhead/Snakefile

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
        counts = "analysis/STAR/STAR_Gene_Counts.csv"
    output:
        limma = "analysis/diffexp/{comparison}/{comparison}.limma.csv",
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.csv",
        deseqSum = "analysis/diffexp/{comparison}/{comparison}.deseq.sum.csv",
        #annotations
        limma_annot = "analysis/diffexp/{comparison}/{comparison}.limma.annot.csv",
        deseq_annot = "analysis/diffexp/{comparison}/{comparison}.deseq.annot.csv",
    params:
        s1=lambda wildcards: ",".join(get_comparison(wildcards.comparison, 1)),
        s2=lambda wildcards: ",".join(get_comparison(wildcards.comparison, 2)),
        gene_annotation = config['gene_annotation']
    message: "Running differential expression analysis using limma and deseq for {wildcards.comparison}"
#    script:
#        "scripts/DEseq.R"

    run:
        shell("Rscript viper/scripts/DEseq.R \"{input.counts}\" \"{params.s1}\" \"{params.s2}\" {output.limma} {output.deseq} {output.limma_annot} {output.deseq_annot} {output.deseqSum} {params.gene_annotation}")

rule fetch_DE_gene_list:
    input:
        deseq_file_list=expand("analysis/diffexp/{comparison}/{comparison}.deseq.csv",comparison=comparisons)
    output:
        csv="analysis/diffexp/de_summary.csv",
        png="analysis/diffexp/de_summary.png"
    message: "Creating Differential Expression summary"
    run:
        deseq_file_string = ' -f '.join(input.deseq_file_list)
        shell("perl viper/scripts/get_de_summary_table.pl -f {deseq_file_string} 1>{output.csv}")
        shell("Rscript viper/scripts/de_summary.R {output.csv} {output.png}")

#Generate volcano plots for each comparison
rule volcano_plot:
    input:
        deseq = "analysis/diffexp/{comparison}/{comparison}.deseq.csv",
    output:
        plot = "analysis/diffexp/{comparison}/{comparison}_volcano.pdf",
        png = "analysis/plots/images/{comparison}_volcano.png"
    message: "Creating volcano plots for Differential Expressions for {wildcards.comparison}"
    run:
        shell("Rscript viper/scripts/volcano_plot.R {input.deseq} {output.plot} {output.png}")

#call snps from the samples
#NOTE: lots of duplicated code below!--ONE SET for chr6 (default) and another
#for genome-wide
#------------------------------------------------------------------------------
# snp calling for chr6 (default)
#------------------------------------------------------------------------------
rule call_snps_chr6:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        ref_fa=config["ref_fasta"],
    output:
        "analysis/snp/{sample}/{sample}.snp.chr6.txt"
    params:
        varscan_jar_path=config["varscan_jar_path"]
    message: "Running varscan for snp analysis for ch6 fingerprint region"
    shell:
        "samtools mpileup -r \"chr6\" -f {input.ref_fa} {input.bam} | awk \'$4 != 0\' | "
        "java -jar {params.varscan_jar_path} pileup2snp - --min-coverage 20 --min-reads2 4 > {output}"

#calculate sample snps correlation using all samples
rule sample_snps_corr_chr6:
    input:
        snps = lambda wildcards: expand("analysis/snp/{sample}/{sample}.snp.chr6.txt", sample=ordered_sample_list)
    output:
        "analysis/snp/snp_corr.chr6.txt"
    message: "Running snp correlations for chr6 fingerprint region"
    run:
        snps = " ".join(input.snps)
        shell("{config[python2]} viper/scripts/sampleSNPcorr.py {snps}> {output}")

rule snps_corr_plot_chr6:
    input:
        snp_corr="analysis/snp/snp_corr.chr6.txt",
        annotFile=config['metasheet'],
    output:
        snp_plot_out="analysis/plots/sampleSNPcorr_plot.chr6.png"
    message: "Running snp analysis for chr6 fingerprint region"
    run:
        shell("Rscript viper/scripts/sampleSNPcorr_plot.R {input.snp_corr} {input.annotFile} {output.snp_plot_out}")

#------------------------------------------------------------------------------
# snp calling GENOME wide (hidden config.yaml flag- 'snp_scan_genome:True'
#------------------------------------------------------------------------------

rule call_snps_genome:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        ref_fa=config["ref_fasta"],
    output:
        "analysis/snp/{sample}/{sample}.snp.genome.txt"
    params:
        varscan_jar_path=config["varscan_jar_path"]
    message: "Running varscan for snp analysis genome wide"
    shell:
        "samtools mpileup -f {input.ref_fa} {input.bam} | awk \'$4 != 0\' | "
        "java -jar {params.varscan_jar_path} pileup2snp - --min-coverage 20 --min-reads2 4 > {output}"

rule sample_snps_corr_genome:
    input:
        snps = lambda wildcards: expand("analysis/snp/{sample}/{sample}.snp.genome.txt", sample=ordered_sample_list)
    output:
        "analysis/snp/snp_corr.genome.txt"
    message: "Running snp analysis genome wide"
    run:
        snps = " ".join(input.snps)
        shell("{config[python2]} viper/scripts/sampleSNPcorr.py {snps}> {output}")

rule snps_corr_plot_genome:
    input:
        snp_corr="analysis/snp/snp_corr.genome.txt",
        annotFile=config['metasheet'],
    output:
        snp_plot_out="analysis/plots/sampleSNPcorr_plot.genome.png"
    message: "Creating snp plot genome wide"
    run:
        shell("Rscript viper/scripts/sampleSNPcorr_plot.R {input.snp_corr} {input.annotFile} {output.snp_plot_out}")

## Perform Correlation analysis between limma diff files
#rule correlation_plot:
#    input:
#        diffiles = expand("analysis/diffexp/{comparison}/{comparison}.deseq.csv", comparison=comparisons),
#        meta = config["metasheet"]
#    output:
#        correlation_plot = "analysis/plots/correlation_plot.pdf",
#        correlation_table = "analysis/plots/correlation_table.csv",
#        upvenn_plot = "analysis/plots/upvenn_plot.pdf",
#        downvenn_plot = "analysis/plots/downvenn_plot.pdf"
#    params:
#        SFnumgenes = config["SFnumgenes"]
#    script:
#        "scripts/correlation_plot.R"


