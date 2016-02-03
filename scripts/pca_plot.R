#Script to generate a PCA (principle components) plot
#Input:
#   rawRPKM_file: a comma-sep file of Gene,Sample1(RPKM),Sample2, ...,SampleN
#   annotation: annotation file in SOME TBD format
#OUTPUT:
#   plot_out: some filename ending in '.pdf' to save the plots
#BASED on garber_analysis.R by Henry Long

# load required packages
library("gplots")
library("ComplexHeatmap")
library("circlize")
library("dendextend")
library("viridis")
library('dplyr')
source('snakemake/scripts/supp_fns.R')

#enable stack trace
#LEN:
options(error = function() traceback(2))

pca_plot <- function(rpkmTable, annotation, plot_out, png_dir) {
    #CONSTANTS
    RPKM_THRESHOLD <- 2.0
    MIN_NUM_SAMPLES_EXPRESSSING_AT_THRESHOLD <- 4
    NUM_GENES_TO_CLUSTER <- 250

    #readin and process newdata
    #newdata <- read.table(args[1], header=T, row.names=1, sep=",")
    newdata <- rpkmTable
    
    #remove MIR RNAs and SNO RNAs    
    #LEN: newdata <- newdata[-grep("MIR[[:digit:]]*",rownames(newdata)), ]
    #LEN: newdata <- newdata[-grep("SNO.*",rownames(newdata)), ]

    #remove genes with no RPKM values or
    #genes where not enough samples meet a minimum threshold
    #LEN: newdata <- (newdata[which(apply(newdata,1,mean)!=0),])
    newdata<-newdata[apply(newdata, 1, function(x) length(x[x>=RPKM_THRESHOLD])>MIN_NUM_SAMPLES_EXPRESSSING_AT_THRESHOLD),]

    #log transform of data
    newdata <- log2(newdata+1)

    #Calculate CVs for all genes (rows)
    mean_rpkm_nolym <- apply(newdata,1,mean)
    var_rpkm_nolym <- apply(newdata,1,var)
    cv_rpkm_nolym <- abs(var_rpkm_nolym/mean_rpkm_nolym)

    #Select out the most highly variable genes into the dataframe 'Exp_data'
    Exp_data <- newdata[order(cv_rpkm_nolym,decreasing=T)[1:NUM_GENES_TO_CLUSTER],]

    #SAVE plot
    pdf(file = plot_out)
    png_counter <- 1
    
    #Standard PCA analysis using all possible annotations
    for(c in colnames(annotation)) {
        ann <- as.matrix(annotation[, c])
        if(length(sort(unique(na.omit(as.vector(ann))))) <7) {
            ClassColors <- cmap(ann)
            
            myColors = ClassColors[ann]
            myColors[which(is.na(myColors))] <- "black"
            png(file=paste(png_dir,"/pca_plot_",png_counter,".png",sep=""))
            png_counter <- png_counter + 1
            pca_output <- make_pca_plots(t(Exp_data), threeD = FALSE, ClassColorings = myColors, pca_title = c, legend_title =  c)
	    dev.off()
	    pca_output <- make_pca_plots(t(Exp_data), threeD = FALSE, ClassColorings = myColors, pca_title = c, legend_title =  c)
		
        }
    }

    #GET percent variances
    pc_var <- signif(100.0 * summary(pca_output)[[6]][2,], digits = 3)
    #scree plot
    png(file=paste(png_dir,"/pca_plot_",png_counter,".png",sep=""))
    barplot(pc_var, ylim=c(0,100),ylab="% variance")
    dev.off()
    barplot(pc_var, ylim=c(0,100),ylab="% variance")
    #Short Summary
    #summary(pca_output)[[6]][,1:3]

    #SAVE graphics
    dev.off()
}

args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
pca_plot_out=args[3]
png_dir=args[4]

#process RPKM file
# Mahesh adding check.names=F so that if there is any - or _ characters, they won't be turned to default '.'
rpkmTable <- read.table(rpkmFile, header=T, check.names=F, row.names=1, sep=",", stringsAsFactors=FALSE, dec=".")
#CONVERT to numeric!
for (n in names(rpkmTable)) {
    rpkmTable[n] <- apply(rpkmTable[n], 1, as.numeric)
}

#PROCESS ANNOTATIONS
tmp_ann <- read.delim(annotFile, sep=",", stringsAsFactors=FALSE)
#REMOVE comp_ columns
#previous Len's code was returning 0 columns if metasheet doesn't contain 'comp_' column 
tmp_ann <- tmp_ann[ , !grepl('comp_*', names(tmp_ann))]
rownames(tmp_ann) <- tmp_ann$SampleName
samples <- intersect(colnames(rpkmTable), rownames(tmp_ann))
tmp_ann <- tmp_ann[samples,-1]
pca_plot(rpkmTable, tmp_ann, pca_plot_out, png_dir)
