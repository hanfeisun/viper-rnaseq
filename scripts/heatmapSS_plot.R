#Script to generate a Sample-Sample heatmaps
#Input:
#   rawRPKM_file: a comma-sep file of Gene,Sample1(RPKM),Sample2, ...,SampleN
#   annotation: annotation file in SOME TBD format
#Output:
#   plot_out: some filename ending in '.pdf' to save the plots
#   ssSpear_out: some filename ending in 'txt' to save the sample-sample
#                   spearman correlations
#BASED on makeclustering_all_samples_HWL.R by Henry Long

# load required packages
library("gplots")
library("ComplexHeatmap")
library("circlize")
library("dendextend")
library("viridis")
library('dplyr')
source('snakemake/scripts/supp_fns.R')

#enable stack trace
options(error = function() traceback(2))

heatmapSS_plot <- function(rpkmTable, annotation, plot_out, ssSpear_out) {
    #CONSTANTS
    RPKM_THRESHOLD <- 2.0
    MIN_NUM_SAMPLES_EXPRESSSING_AT_THRESHOLD <- 4
    NUM_GENES_TO_CLUSTER <- 250

    #readin and process newdata
    #newdata <- read.table(rawRPKM_file, header=T, row.names=1, sep=",")
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

    #Calc. spearman correlation
    cordata <- cor(Exp_data, method="spearman")

    #save the spearman correlation as txt; save plots
    pdf(file = plot_out)
    #data.matrix.filename <- "XXX_foo.txt"
    write.table(cordata, file=ssSpear_out, quote=F, col.names = NA, sep="\t")

    #make SS (sample-sample) heatmap
    ma_nolym <- max(cordata)
    mi_nolym <- min(cordata)
    my.breaks_nolym<-c(mi_nolym,seq(mi_nolym + 0.01, ma_nolym - 0.01,length.out=99),ma_nolym)

    test<-heatmap.2(as.matrix(cordata), Colv = "Rowv", revC = T,
                    distfun = function(x) dist(x,method = 'euclidean'),
                    hclustfun = function(x) hclust(x,method = 'ward.D2'),
                    breaks=my.breaks_nolym,trace="none",scale="none",
                    col=bluered(100),labCol=F, cexRow=0.75,
                    key=TRUE, margins=c(2,6),
                    #main = paste(Project_Name,"Sample-Sample Correlation"),
                    main = "Sample-Sample Correlation",
                    #xlab = paste(NUM_GENES_TO_CLUSTER, " genes ",Sys.Date())
                    xlab = paste(NUM_GENES_TO_CLUSTER, " genes ")
                    )

    #LEN: HYP- This is generating matrix_1
    ss_col = colorRamp2(seq(min(cordata), max(cordata), length = 3), c("blue", "#EEEEEE", "red"))
    ht_list <- Heatmap(cordata, name="sprmanCorr", col = ss_col)
    ht_list <- make_complexHeatmap_annotation(ht_list, annotation)
    draw(ht_list)
    
    #SAVE graphics
    dev.off()
}

args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
ss_plot_out=args[3]
ss_txt_out=args[4]

#process RPKM file
rpkmTable <- read.table(rpkmFile, header=T, row.names=1, sep=",", stringsAsFactors=FALSE, dec=".")
for (n in names(rpkmTable)) {
    #CONVERT to numeric!
    rpkmTable[n] <- apply(rpkmTable[n], 1, as.numeric)
    #replace NA
    #na.omit(rpkmTable)
}

#PROCESS ANNOTATIONS
tmp_ann <- read.delim(annotFile, sep=",", stringsAsFactors=FALSE)
#REMOVE comp_ columns
tmp_ann <- tmp_ann[ , -grep('comp_*', names(tmp_ann))]
rownames(tmp_ann) <- tmp_ann$SampleName
samples <- intersect(colnames(rpkmTable), rownames(tmp_ann))
tmp_ann <- tmp_ann[samples,-1]
#print(str(tmp_ann))
heatmapSS_plot(rpkmTable, tmp_ann, ss_plot_out, ss_txt_out)
