suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
### NEED MOUSE DB
suppressMessages(library("GOstats"))
suppressMessages(library("edgeR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggalt"))
suppressMessages(library("scales"))


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

    ##### LOLLIPOP, AWAITING ggalt 0.3.0.0

    #define limits of annotation label - note that since the function 'strwrap' tries to make breaks at word boundaries it is possible for the label to wrap over (go.term.height + 1) lines
    #go.term.width = 60
    #go.term.height = 2
    #max.length = go.term.height*go.term.width

    #trim white space and then determine which GO terms in the summary are long enough to be folded
    #df$Term <- trimws(df$Term)
    #to_fold <- which(lapply(df$Term, nchar) > go.term.width)

    #df$Term[to_fold] <- sapply(df$Term[to_fold],function(x) ifelse(nchar(x)>max.length, substr(x,1,max.length),x) )
    #df$Term[to_fold] <- sapply(df$Term[to_fold],function(x) paste(strwrap(x,width = go.term.width),collapse = "\n"))

    # Use ggplt2 to make a custom 'lollipop' chart of the results
    # Arial Narrow is a good font to put as many characters as possible into a small space
    # Exact font sizes and faces can be tweaked futher
    #pdf("/mnt/cfce-stor1/home/mgc31/code/viperproject/test.pdf")

    #gg <- ggplot(df, aes(y=reorder(Term, logpval), x=logpval))
    #gg <- gg + geom_lollipop(point.colour="steelblue", point.size=5, horizontal=TRUE)
    #gg <- gg + labs(x="- Log(P-value)", y=NULL, title="GO Term analysis")
    #gg <- gg + theme_minimal(base_family="Arial Narrow")
    #gg <- gg + theme(panel.grid.major.y=element_blank())
    #gg <- gg + theme(panel.grid.minor=element_blank())
    #gg <- gg + theme(axis.text.y=element_text(debug = FALSE, size=14, lineheight = 0.9, margin=margin(0,-20,0,20)))
    #gg <- gg + theme(axis.text.x=element_text(debug = FALSE, size=16, margin=margin(-10,0,20,0)))
    #gg <- gg + theme(axis.title.x=element_text(debug = FALSE, size=18, margin=margin(0,0,20,0)))
    #gg <- gg + theme(plot.title = element_text(size=28, margin = margin(10, 0, 20, 0)))
    #gg <- gg + geom_vline(xintercept = 0)

    #dev.off()

    

    ## Create title for plot
    temptitle = tail(unlist(strsplit(goterm_pdf, split="/")), n=1)
    temptitle = head(unlist(strsplit(temptitle, split="[.]")), n=1)
    title = paste(temptitle, "_Top_", numgoterms, "_GOterms", sep="")

    go_bar_plot <- ggplot(df[numgoterms:1,], aes(factor(Term, levels=unique(Term)), logpval)) + 
      ylab("-log(Pvalue)") + xlab("Go term") +
      geom_bar(stat = "identity", fill="palegreen3") + theme_bw(base_size = 12) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40, indent = 2),"\n") + coord_flip() + ggtitle(title)

    ggsave(goterm_pdf, width=11, height=8.5, unit="in")
    ggsave(goterm_png, width=10, height=8, unit="in")

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
