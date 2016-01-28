#libraries

volcano_plot_f <- function(deseq_results, plot_out) {
    #CREATE pdf output file
    pdf(file = plot_out)

    #FROM: http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
    res <- read.table(deseq_results, header=TRUE)

    #get comparisonName
    comparisonName = strsplit(deseq_results, "/")[[1]][3]

    # Make a basic volcano plot
    #LEN NOTE: the original graph was cutting things off, in terms of log2FC
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=comparisonName, col=rgb(1,1,1,0)))

    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    #significant = blue, non-sig = red, 30% alpha
    with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col=rgb(70,131,180,75, maxColorValue=255))) #col="blue"))
    with(subset(res, padj>=.05), points(log2FoldChange, -log10(pvalue), pch=20, col=rgb(255,0,0,75,  maxColorValue=255)))
    #with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pval), pch=20, col=rgb(0,0,0,127.5,  maxColorValue=255)))#col="orange"))
    #with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pval), pch=20, col=rgb(0,255,0,127.5,  maxColorValue=255)))#"green"))


    # Label points with the textxy function from the calibrate plot
    library(calibrate)
    #LEN NOTE: this cutoff is too liberal
    #with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pval), labs=id, cex=.8))
    #BETTER
    #with(subset(res, padj<10e-7), textxy(log2FoldChange, -log10(pval), labs=id, cex=.8))
    
    #LEN NOTE: would be best if could label top 100 sig genes.
    topSig <- res[order(res$pvalue),]
    #print(head(topSig))
    with(head(topSig, 100), textxy(log2FoldChange, -log10(pvalue), labs=id, cex=.5))

    dev.off()
}

args <- commandArgs( trailingOnly = TRUE )
deseq = args[1]
out = args[2]

volcano_plot_f(deseq, out)
