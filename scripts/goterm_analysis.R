suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
### NEED MOUSE DB
suppressMessages(library("GOstats"))
suppressMessages(library("edgeR"))

goterm_analysis_f <- function(deseq_file, goterm_csv,goterm_pdf,goterm_png) {

    ### PARAMS:
    adjpvalcutoff = 0.05 ## for genes that will go into gostats for goterm analysis
    numgoterms = 20 ## number of go terms to bar chart
    reference = "hg19"

    ## Read in detable
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]
    
    ## Append ENSEMBL and ENTREZ IDs from loaded in db
    if (reference == "hg19") {IDdb = org.Hs.eg.db}
    #if (reference == "mm9") {IDdb = ##########}

    detable$ensembl <- mapIds(IDdb,
                         keys=rownames(detable),
                         column="ENSEMBL",
                         keytype="SYMBOL",
                         multiVals="first")
    detable$entrez <- mapIds(IDdb,
                         keys=rownames(detable),
                         column="ENTREZID",
                         keytype="SYMBOL",
                         multiVals="first")

    adjpvalcutoff = 0.05

    ## Select genes that pass the adjPval cutoff and select those entrez IDs as pop, set rest as universe.
    topgenes <- subset(detable, detable$padj < adjpvalcutoff)
    selectedIDs = topgenes$entrez
    universeIDs = detable$entrez

    ## Run GOstats
    goParams <- new("GOHyperGParams",
                    geneIds = selectedIDs,
                    universeGeneIds = universeIDs,
                    annotation = "org.Hs.eg.db",
                    ontology = "BP", ## CC, BP, MF
                    pvalueCutoff = adjpvalcutoff,
                    conditional = TRUE,
                    testDirection = "over")

    goResults <- hyperGTest(goParams)

    ## Summary table has columns: GOBPID, Pvalue, OddsRatio, ExpCount, Count, Size, Term.
    df = summary(goResults)
    #df$percent = df$Count/length(selectedIDs)
    df$logpval = -log(df$Pvalue)

    ## Write out Results
    write.table(df, file = goterm_csv, col.names=T, row.names=F, quote=F, sep=",")
    
    numgoterms = 20

    # Core wrapping functions for barplot mapping
    wrap.it <- function(x, len) {sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)}
    wrap.labels <- function(x, len) {if (is.list(x)) {lapply(x, wrap.it, len)} else {wrap.it(x, len)}}

    ## Function for rounding axes to upper 5 limit
    mround <- function(x,base) {base*ceiling(x/base)}

    ## Make GOterm table labels
    wr.lap = wrap.labels(df$Term[numgoterms:1], 35)
    
    ## Create title for plot
    temptitle = tail(unlist(strsplit(goterm_pdf, split="/")), n=1)
    temptitle = head(unlist(strsplit(temptitle, split="[.]")), n=1)
    title = paste(temptitle, "_Top_", numgoterms, "_GOterms", sep="")

    ## Plot out in pdf
    pdf(goterm_pdf, width=11,height=8.5)
    par(mar=c(5,10,2,1))
    
    barplot(df$logpval[numgoterms:1], names.arg = wr.lap,
            width = 3, space = 0.2,
            main = title, horiz = TRUE, xlab = "-log(pvalue)",
            xlim = c(0,mround(max(df$logpval[numgoterms:1]),5)), xpd = TRUE, las = 2, cex.names = 0.7)

    dev.off()

    ## Plot out in png
    png(goterm_png, width = 10, height = 8, unit="in",res=300)
    par(mar=c(5,10,2,1))

    barplot(df$logpval[numgoterms:1], names.arg = wr.lap,
            width = 3, space = 0.2,
            main = title, horiz = TRUE, xlab = "-log(pvalue)",
            xlim = c(0,mround(max(df$logpval[numgoterms:1]),5)), xpd = TRUE, las = 2, cex.names = 0.7)

    dev.off()

}

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
goterm_csv = args[2]
goterm_pdf = args[3]
goterm_png = args[4]

goterm_analysis_f(deseq_file, goterm_csv,goterm_pdf,goterm_png)



#goterm_anlysis_f(
#    snakemake@input[[""]],
#    snakemake@params[[""]],
#    snakemake@output[[""]]
#    )
