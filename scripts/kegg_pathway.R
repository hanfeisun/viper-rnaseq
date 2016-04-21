kegg_pathway_f<- function(deseq_file) {

        detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)


    }

args <- commandArgs( trailingOnly = TRUE )
deseq <- file= args[1]
pdf <- file = args[2]
png <- file = args[3]

deseq_file = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/D538G_WM_NoDOXvDOX/D538G_WM_NoDOXvDOX.deseq.annot.csv"


kegg_pathway_f(deseq_file)






#kegg_pathway_f(
#    snakemake@input[[""]],
#    snakemake@params[[""]],
#    snakemake@output[[""]]
#    )
