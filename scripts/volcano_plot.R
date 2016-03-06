#libraries
library(calibrate)

volcano_plot_f <- function(deseq_results, pdf_file, png_file, makePDF = TRUE) {
    #NOTE: if makePDF is FALSE, generate PNG file

    if (makePDF) {
        #CREATE pdf as output file
        pdf(file = pdf_file)
    } else {
        #CREATE png as output file
        png(file=png_file, width = 8, height = 8, unit="in",res=300)
    }

    #FROM: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
    res <- read.table(deseq_results, header=TRUE, sep=",")

    #get comparisonName
    comparisonName = strsplit(deseq_results, "/")[[1]][3]

    # Make a basic volcano plot
    #LEN NOTE: the original graph was cutting things off, in terms of log2FC
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=comparisonName, col="gray"))

    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    #significant = blue, non-sig = red, 30% alpha
    with(subset(res, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col=rgb(70,131,180,75, maxColorValue=255))) #col="blue"))
    with(subset(res, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col=rgb(255,0,0,75,  maxColorValue=255))) #RED

    # Label points with the textxy function from the calibrate plot
    #LEN NOTE: this cutoff is too liberal--I re-enabled it, and it's crowded!
    with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=id, cex=.4))
    #BETTER
    #with(subset(res, padj<10e-7), textxy(log2FoldChange, -log10(pval), labs=id, cex=.8))

    ##LABEL the top 100 gene (by pval)
    #topSig <- res[order(res$pvalue),]
    ##print(head(topSig))
    #with(head(topSig, 100), textxy(log2FoldChange, -log10(pvalue), labs=id, cex=.5))
    
    dev.off()

    #--------------------------------------------------------------------------
    #REMOVED: DUPLICATED CODE
    #--------------------------------------------------------------------------
}

args <- commandArgs( trailingOnly = TRUE )
deseq = args[1]
pdf_file = args[2]
png_file = args[3]
volcano_plot_f(deseq, pdf_file, png_file, makePDF = TRUE) #MAKE PDF
volcano_plot_f(deseq, pdf_file, png_file, makePDF = FALSE) #MAKE PNG

