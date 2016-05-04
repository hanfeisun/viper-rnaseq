suppressMessages(library("dplyr"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("gage"))
suppressMessages(library("gageData"))
suppressMessages(library("pathview"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("XML"))
data(kegg.sets.hs)
data(sigmet.idx.hs)

## "cleans up the gene sets
#kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

#deseq_file = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/MCF7_PvTAMR/MCF7_PvTAMR.deseq.csv"
#reference = "hg19"

#kegg_table = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/MCF7_PvTAMR/MCF7_PvTAMR.kegg.txt"
#kegg_dir= "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/MCF7_PvTAMR/kegg_pathways/"
#gsea_table = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/MCF7_PvTAMR/MCF7_PvTAMR.gsea.txt"
#gsea_pdf = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/MCF7_PvTAMR/MCF7_PvTAMR.gsea.pdf"
#temp_dir = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/diffexp/MCF7_PvTAMR/temp/"

kegg_pathway_f<- function(deseq_file, kegg_dir,reference,temp_dir, kegg_table,gsea_table,gsea_pdf) {

    mainDir = substr(kegg_dir, 1, nchar(kegg_dir)-14)
    dir.create(file.path(mainDir, "kegg_pathways/"), showWarnings = FALSE)

    ## Read in deseq table
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]

    #detable = subset(detable, abs(detable[,3]) > 1 & detable[,7] < 0.05)

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

    ## Setting up gage input, needs the log2fc with the entrez id
    gageinput = detable$log2FoldChange
    names(gageinput) = detable$entrez

    ## Run gage
    keggres = gage(gageinput, gsets = kegg.sets.hs, same.dir=TRUE)
    
    kegg_output = keggres$greater
    kegg_output = cbind(rownames(kegg_output), kegg_output)
    colnames(kegg_output)[1] = "Kegg_pathway"
    write.table(kegg_output, file = kegg_table, quote=F, col.names=TRUE, row.names=FALSE, sep="\t")
    
    ## Get the pathways
    numpathways = 5

    keggrespathways = keggres$stats
    keggrespathways = keggrespathways[order(-abs(keggrespathways[,1])),]
    keggrespathways = rownames(keggrespathways)[1:numpathways]
    
    #keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
    #  tbl_df() %>%
    #  filter(row_number()<=numpathways) %>%
    #  .$id %>%
    #  as.character()

    keggresids = substr(keggrespathways, start=1, stop=8)

    ## Plot using pathview
    #plot_pathway = function(pid) pathview(gene.data=gageinput, pathway.id=pid, species="hsa", new.signature=FALSE, kegg.dir=temp_dir)
    #tmp = sapply(keggresids, function(pid) pathview(gene.data=gageinput, pathway.id=pid, species="hsa", kegg.dir=temp_dir))

    normwd = getwd()
    setwd(temp_dir)
    
    for ( i in 1:numpathways) {
        pvout <- pathview(gene.data=gageinput,              ## Gene list
                          pathway.id=keggresids[i],         ## Which pathway
                          species = "hsa",                  ## Species
                          #limit = list(gene=max(abs(gageinput)),cpd=1),
                          #kegg.dir = temp_dir               ## Save directory
                          low = list(gene = "blue", cpd = "purple"),
                          mid = list(gene = "gray", spd = "gray"),
                          high = list(gene = "red", cpd = "gray")  ## Color scale
                     )
    }
    setwd(normwd)

    ## Renaming files
    # Create variable with keggrespathways sorted and pull out name
    sortkeggrespathways = sort(keggrespathways)
    newnames = substr(sortkeggrespathways, 10, nchar(sortkeggrespathways))
    newnames = gsub(" ", "_", newnames)

    # Read in the list of made png files
    png_files <- list.files(temp_dir, pattern=glob2rx("*.pathview.png"))

    file.rename(paste0(temp_dir,png_files), paste0(kegg_dir, newnames, ".png"))

    # Repeat for xml files
    xml_files <- list.files(temp_dir, pattern=glob2rx("*.xml"))

    file.rename(paste0(temp_dir,xml_files), paste0(kegg_dir, newnames, ".xml"))

    ## GSEA Analysis
    
    #upgenes = subset(gageinput, gageinput > 0)
    #upgenes = sort(upgenes, decreasing = TRUE)
    #downgenes = subset(gageinput, gageinput < 0)
    #downgenes = -sort(downgenes)

    gseainput = sort(gageinput, decreasing=TRUE)

    fullgsea <- gseKEGG(geneList = gseainput,
                        organism     = "human",
                        nPerm        = 100,
                        minGSSize    = 1,
                        pvalueCutoff = 0.99,
                        verbose      = FALSE,
                        use_internal_data = FALSE)


    gsea_data = summary(fullgsea)
    gsea_data = gsea_data[order(-abs(gsea_data$NES)),]
    write.table(gsea_data, file = gsea_table, quote=FALSE, sep= "\t", row.names=FALSE, col.names=TRUE)

    pdf(gsea_pdf)
    plot.new()
    mtext("GSEA_plots")
    for ( i in 1:10 ) {
        gseaplot(fullgsea, geneSetID = gsea_data[i,1])
        mtext(gsea_data[i,2])
    }
    dev.off()
    
    }

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
kegg_dir = args[2]
reference = args[3]
temp_dir = args[4]
kegg_table = args[5]
gsea_table = args[6]
gsea_pdf = args[7]


kegg_pathway_f(deseq_file, kegg_dir,reference,temp_dir, kegg_table,gsea_table,gsea_pdf)






#kegg_pathway_f(
#    snakemake@input[[""]],
#    snakemake@params[[""]],
#    snakemake@output[[""]]
#    )
