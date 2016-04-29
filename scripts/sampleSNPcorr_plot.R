#Script to generate a Sample-Sample snp correlation heatmaps
#Input:
#Output:
#BASED on makeclustering_all_samples_HWL.R by Henry Long

# load required packages
library("gplots")
library("ComplexHeatmap")
library("circlize")
library("dendextend")
library("viridis")
library('dplyr')
source('viper/scripts/supp_fns.R')

#enable stack trace
options(error = function() traceback(2))

#NOTE: this fn is called twice, once to generate png and another to make pdf
snp_corr_plot <- function(snpCorrMatrix, annotation, plot_out, isPNG) {
    cordata <- snpCorrMatrix

    ## We want to only work with the samples that are in the meta file, so we are only selecting the count columns that are in the meta file
    cordata <- cordata[, rownames(annotation)]
    cordata <- cordata[rownames(annotation),]
    
    #save the spearman correlation as txt; save plots
    if (isPNG) {
        png(file = plot_out)
    } else {
        pdf(file = plot_out)
    }

    #LEN: Henry's updated plotting code
    dataset = cordata
    breakpoint = 0.75
    my_breaks = c(seq(0,breakpoint,length=20),seq(breakpoint,1.0,length=180))
    my_palette <- colorRampPalette(c("white", "red"))(n = 199)

    heatmap.2(as.matrix(dataset),
                     dendrogram="none",
                     trace = 'none',
                     breaks = my_breaks,
                     col = my_palette,
                     key = FALSE,
                     # block sepration
                     colsep = 1:ncol(dataset),
                     rowsep = 1:nrow(dataset),
                     sepcolor="white",
                     sepwidth=c(0.02,0.02),
                     margins=c(10,10),
                     main = "SNP Correlation Plot")
    #draw(graph)
    dev.off()

    #LEN: OLD plotting code
    
    ## #make SS (sample-sample) heatmap
    ## ma_nolym <- max(cordata)
    ## mi_nolym <- min(cordata)
    ## my.breaks_nolym<-c(mi_nolym,seq(mi_nolym + 0.01, ma_nolym - 0.01,length.out=99),ma_nolym)
    
    ## #ha1 <- make_complexHeatmap_annotation(annotation)
    ## graph2 <-Heatmap(t(as.matrix(cordata)), name="scale",
    ##                  col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
    ##                  #column_dend_height = unit(2, "cm"),
    ##                  #heatmap_legend_param = list(title = "exp. level"),
    ##                  column_title = "Sample-Sample Correlation",
    ##                  #row_title = "Samples",
    ##                  show_row_names = TRUE, show_column_names = TRUE,
    ##                  #row_names_max_width = unit(3, "mm"),
    ##                  row_names_gp = gpar(fontsize = 12),
    ##                  column_names_gp = gpar(fontsize = 12),
    ##                  #SETTING the diagonal order
    ##                  cluster_rows = TRUE, cluster_columns=TRUE,
    ##                  row_order=1:nrow(cordata), #coloumn_order=1:ncol(cordata),
    ##                  show_column_dend=FALSE, show_row_dend=FALSE,
    ##                  clustering_method_rows="ward.D2",
    ##                  clustering_method_columns="ward.D2",
    ##                  clustering_distance_rows="euclidean",
    ##                  clustering_distance_columns="euclidean",
    ##                  show_heatmap_legend = TRUE,
    ##                  #row_dend_width = unit(5, "mm"),
    ##                  #width=unit(60,"cm"),
    ##                  #top_annotation=ha1,
    ##                  )
    ##draw(graph2)
    ##dev.off()
    
}

args <- commandArgs( trailingOnly = TRUE )
snpCorrFile=args[1]
annotFile=args[2]
snp_corr_plot_out=args[3]
snp_corr_plot_pdf=args[4]

#READ in corr. file
snpCorrMat <- read.table(snpCorrFile, header=TRUE, sep="\t", row.names=1)

#NOTE: in the snpCorrMatrix, sample i.e. column and row names are in the
#form SAMPLEXXX.snp.chr6 or SAMPLEXXX.snp.[something]
#WE need to EXTRACT out SAMPLEXXX from this string -> sampleNames
sampleNames <- sapply(colnames(snpCorrMat),
                      function(x) substring(x, 1, regexpr('.snp',x) - 1))
colnames(snpCorrMat) <- sampleNames
rownames(snpCorrMat) <- sampleNames

#PROCESS ANNOTATIONS
tmp_ann <- read.delim(annotFile, sep=",", stringsAsFactors=FALSE)
#REMOVE comp_ columns
tmp_ann <- tmp_ann[ , !grepl('comp_*', names(tmp_ann))]

#convert numerical annotations to numbers/floats
for (col in colnames(tmp_ann)) {
    #IS it a valid number?--test first value in col
    if(attr(regexpr("^\\-?\\d+\\.\\d+$",tmp_ann[1,col]), "match.length") > 0){
        #print(apply(as.matrix(tmp_ann[,col]), 2, as.numeric))
        tmp_ann[,col] <- as.vector(apply(as.matrix(tmp_ann[,col]), 2, as.numeric))
    }
}

rownames(tmp_ann) <- tmp_ann[,1]
samples <- intersect(colnames(snpCorrMat), rownames(tmp_ann))
tmp_ann <- tmp_ann[samples,-1]
#print(str(tmp_ann))

#GENERATE png
snp_corr_plot(snpCorrMat, tmp_ann, snp_corr_plot_out, T)
#GENERATE png
snp_corr_plot(snpCorrMat, tmp_ann, snp_corr_plot_pdf, F)
