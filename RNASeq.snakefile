# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
from collections import defaultdict

configfile: "/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/snakemake/mm9.yaml"

workdir: config["directory"]

strand_command=""
rRNA_strand_command=""

if( config["stranded"] ):
    strand_command="--outFilterIntronMotifs RemoveNoncanonical"
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
        file_info[info[1]].append(info[0])
        if info[1] not in ordered_sample_list:
            ordered_sample_list.append(info[1])
        if ( not run_fusion and "_R2_" in info[0] ):
            run_fusion=True

if( run_fusion ):
    if( config["stranded"] ):
        strand_command = " --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000"
    else:
        strand_command = " --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMstrandField intronMotif"


def get_fastq(wildcards):
    return [os.path.join("concat_per_sample_fastq", f) for f in file_info[wildcards.sample]]


def fusion_output(wildcards):
    fusion_out_files = []
    if run_fusion:
        for sample in file_info.keys():
            fusion_out_files.append( "STAR_Fusion/" + sample + "/" + sample + ".fusion_candidates.final" )
    return fusion_out_files

def insert_size_output(wildcards):
    insert_size_out_files = []
    if run_fusion:
        for sample in file_info.keys():
            insert_size_out_files.append( "RSeQC/insert_size/" + sample + "/" + sample + ".histogram.pdf" )
    return insert_size_out_files

def rRNA_metrics(wildcards):
    if config["star_rRNA_index"] is not None:
        return "STAR_rRNA/STAR_rRNA_Align_Report.csv"

rule target:
    input:
        expand( "cufflinks/{K}/{K}.genes.fpkm_tracking", K=ordered_sample_list ),
        "STAR_snakemake/STAR_Align_Report.csv",
        "STAR_snakemake/STAR_Align_Report.png",
        "STAR_snakemake/STAR_Gene_Counts.csv",
        "cufflinks/Cuff_Gene_Counts.csv",
        expand( "RSeQC/read_distrib/{sample}.txt", sample=ordered_sample_list ),
        "RSeQC/read_distrib/read_distrib.png",
        expand( "RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png", sample=ordered_sample_list ),
        "RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        expand( "RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf", sample=ordered_sample_list ),
        fusion_output,
        insert_size_output,
        rRNA_metrics

rule run_STAR:
    input:
        get_fastq
    output:
        bam="STAR_snakemake/{sample}/{sample}.sorted.bam",
        counts="STAR_snakemake/{sample}/{sample}.counts.tab",
        log_file="STAR_snakemake/{sample}/{sample}.Log.final.out"
    params:
        stranded=strand_command,
        prefix=lambda wildcards: "STAR_snakemake/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_index]}"
        " --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.prefix}."
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000 --quantMode GeneCounts"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"
        " && /usr/bin/samtools index {output.bam}"


rule generate_STAR_report:
    input:
        star_log_files=expand( "STAR_snakemake/{sample}/{sample}.Log.final.out", sample=ordered_sample_list ),
        star_gene_count_files=expand( "STAR_snakemake/{sample}/{sample}.counts.tab", sample=ordered_sample_list )
    output:
        csv="STAR_snakemake/STAR_Align_Report.csv",
        png="STAR_snakemake/STAR_Align_Report.png",
        gene_counts="STAR_snakemake/STAR_Gene_Counts.csv"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/map_stats.R {output.csv} {output.png}" )
        shell( "perl /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/rna_seq/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

rule run_cufflinks:
    input:
        "STAR_snakemake/{sample}/{sample}.sorted.bam"
    output:
        "cufflinks/{sample}/{sample}.genes.fpkm_tracking"
    threads: 4
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/cufflinks-2.2.1.Linux_x86_64/cufflinks -o cufflinks/{wildcards.sample} -p {threads} -G {config[gtf_file]} {input}"
        " && mv cufflinks/{wildcards.sample}/genes.fpkm_tracking {output}"
        " && mv cufflinks/{wildcards.sample}/isoforms.fpkm_tracking cufflinks/{wildcards.sample}/{wildcards.sample}.isoforms.fpkm_tracking"

rule generate_cuff_matrix:
    input:
        cuff_gene_fpkms=expand( "cufflinks/{sample}/{sample}.genes.fpkm_tracking", sample=ordered_sample_list )
    output:
        "cufflinks/Cuff_Gene_Counts.csv"
    run:
        fpkm_files= " -f ".join( input.cuff_gene_fpkms )
        shell( "perl /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/rna_seq/raw_and_fpkm_count_matrix.pl -c -f {fpkm_files} 1>{output}" )


rule run_STAR_fusion:
    input:
        bam="STAR_snakemake/{sample}/{sample}.sorted.bam" #just to make sure STAR output is available before STAR_Fusion
    output:
        "STAR_Fusion/{sample}/{sample}.fusion_candidates.final"
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR-Fusion/STAR-Fusion --chimeric_junction STAR_snakemake/{wildcards.sample}/{wildcards.sample}.Chimeric.out.junction "
        "--chimeric_out_sam STAR_snakemake/{wildcards.sample}/{wildcards.sample}.Chimeric.out.sam --out_prefix STAR_Fusion/{wildcards.sample}/{wildcards.sample}."


rule read_distrib_qc:
    input:
        "STAR_snakemake/{sample}/{sample}.sorted.bam"
    output:
        "RSeQC/read_distrib/{sample}.txt"
    shell:
        "/apps/python-2.7.9/bin/python /zfs/cores/mbcf/mbcf-storage/devel/umv/software/RSeQC/RSeQC-2.6.2/scripts/read_distribution.py"
        " --input-file={input}"
        " --refgene={config[bed_file]} 1>{output}"

rule read_distrib_qc_matrix:
    input:
        read_distrib_files=expand( "RSeQC/read_distrib/{sample}.txt", sample=ordered_sample_list )
    output:
        matrix="RSeQC/read_distrib/read_distrib.matrix.tab",
        png="RSeQC/read_distrib/read_distrib.png"
    run:
        file_list_with_flag = " -f ".join( input.read_distrib_files )
        shell( "perl /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/RSeQC/read_distrib_matrix.pl -f {file_list_with_flag} 1>{output.matrix}" )
        shell( "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/R/R-3.2.2/bin/Rscript /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/RSeQC/read_distrib.R {output.matrix} {output.png}" )


rule down_sample:
    input:
        "STAR_snakemake/{sample}/{sample}.sorted.bam"
    output:
        "RSeQC/gene_body_cvg/downsample/{sample}.downsample.sorted.bam"
    shell:
        "/ifs/rcgroups/mbcf/umv/software/jdk1.8.0_45/bin/java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/picard-tools-1.115/DownsampleSam.jar INPUT={input} OUTPUT={output}"
        " PROBABILITY=0.1"
        " && /usr/local/bin/samtools index {input}"

rule gene_body_cvg_qc:
    input:
        "STAR_snakemake/{sample}/{sample}.sorted.bam"
    output:
        "RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.curves.png",
        "RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r"
    shell:
        "/apps/python-2.7.9/bin/python /zfs/cores/mbcf/mbcf-storage/devel/umv/software/RSeQC/RSeQC-2.6.2/scripts/geneBody_coverage.py -i {input} -r {config[bed_file]}"
        " -f png -o RSeQC/gene_body_cvg/{wildcards.sample}/{wildcards.sample}"


rule plot_gene_body_cvg:
    input:
        expand("RSeQC/gene_body_cvg/{sample}/{sample}.geneBodyCoverage.r", sample=ordered_sample_list )
    output:
        rscript="RSeQC/gene_body_cvg/geneBodyCoverage.r",
        png="RSeQC/gene_body_cvg/geneBodyCoverage.heatMap.png",
        png_curves="RSeQC/gene_body_cvg/geneBodyCoverage.curves.png"
    shell:
        "perl /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/RSeQC/plot_gene_body_cvg.pl --rfile {output.rscript} --png {output.png} --curves_png {output.png_curves} {input}"
        " && /zfs/cores/mbcf/mbcf-storage/devel/umv/software/R/R-3.2.2/bin/Rscript {output.rscript}"

rule junction_saturation:
    input:
        "STAR_snakemake/{sample}/{sample}.sorted.bam"
    output:
        "RSeQC/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf"
    shell:
        "/apps/python-2.7.9/bin/python /zfs/cores/mbcf/mbcf-storage/devel/umv/software/RSeQC/RSeQC-2.6.2/scripts/junction_saturation.py -i {input} -r {config[bed_file]}"
        " -o RSeQC/junction_saturation/{wildcards.sample}/{wildcards.sample}"


rule collect_insert_size:
    input:
        "STAR_snakemake/{sample}/{sample}.sorted.bam"
    output:
        "RSeQC/insert_size/{sample}/{sample}.histogram.pdf"
    shell:
        "/ifs/rcgroups/mbcf/umv/software/jdk1.8.0_45/bin/java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/picard-tools-1.115/CollectInsertSizeMetrics.jar"
        " H={output} I={input} O=RSeQC/insert_size/{wildcards.sample}/{wildcards.sample} R={config[ref_fasta]}"


rule run_rRNA_STAR:
    input:
        get_fastq
    output:
        bam="STAR_rRNA/{sample}/{sample}.sorted.bam",
        log_file="STAR_rRNA/{sample}/{sample}.Log.final.out"
    params:
        stranded=rRNA_strand_command,
        prefix=lambda wildcards: "STAR_rRNA/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample)
    threads: 8
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN {threads} --genomeDir {config[star_rRNA_index]}"
        " --readFilesIn {input} --readFilesCommand zcat --outFileNamePrefix {params.prefix}."
        "  --outSAMmode Full --outSAMattributes All {params.stranded} --outSAMattrRGline {params.readgroup} --outSAMtype BAM SortedByCoordinate"
        "  --limitBAMsortRAM 45000000000"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && /usr/bin/samtools index {output.bam}"


rule generate_rRNA_STAR_report:
    input:
        star_log_files=expand( "STAR_rRNA/{sample}/{sample}.Log.final.out", sample=ordered_sample_list )
    output:
        csv="STAR_rRNA/STAR_rRNA_Align_Report.csv",
        png="STAR_rRNA/STAR_rRNA_Align_Report.png"
    run:
        log_files = " -l ".join( input.star_log_files )
        shell( "perl /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/map_stats.R {output.csv} {output.png}" )
