suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))


kegg_pathway_f<- function(deseq_file) {

    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]

    ## Append ENSEMBL and ENTREZ IDs from loaded in db
    if (reference == "hg19") {IDdb = org.Hs.eg.db}
    #if (reference == "mm9") {IDdb = ##########}
        
    detable$ensembl = mapIds(IDdb,
                        keys=row.names(detable),
                        column="ENSEMBL",
                        keytype="SYMBOL",
                        multiVals="first")
    detable$entrez = mapIds(IDdb,
                        keys=row.names(detable),
                        column="ENTREZID",
                        keytype="SYMBOL",
                        multiVals="first")
    detable$description = mapIds(IDdb,
                        keys=row.names(detable),
                        column="GENENAME",
                        keytype="SYMBOL",
                        multiVals="first")


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
