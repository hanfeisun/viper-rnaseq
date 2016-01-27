#Script to generate a Sample-Feather heatmaps
#Input:
#   rawRPKM_file: a comma-sep file of Gene,Sample1(RPKM),Sample2, ...,SampleN
#   annotation: annotation file in SOME TBD format
#OUTPUT:
#   plot_out: some filename ending in '.pdf' to save the plots
#   sfCorr_out: a .txt file to save sample-feature correlations
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

make_complexHeatmap_annotation <- function(ht_list, annotation){
    MIN_UNIQUE <- 6
    global_gp = gpar(fontsize = 8)
    title_gp = gpar(fontsize = 8, fontface = "bold")
    
    for(c in colnames(annotation)) {
        ann <- as.matrix(annotation[, c])
        if(length(sort(unique(na.omit(as.vector(ann))))) < MIN_UNIQUE) {
            col <- cmap(ann)
        } else {
            #LEN: change this!
            col <- cmap(ann)
            #BELOW causes a bug!--a not a range bug
            #col <- colorRamp2(c(min(ann, na.rm = TRUE), max(ann, na.rm = TRUE)), c("white", "red"))
        }
        jj<- NULL
        #INSTEAD of heatmaps, we are going to add heatmapAnnotations!
        ## foo <- Heatmap(
        ##     t(as.matrix(ann)),
        ##     na_col = "black",
        ##     name = gsub("_", " ", c),
        ##     col = col,
        ##     width = unit(3, "mm"),
        ##     column_names_gp = global_gp,
        ##     column_title_gp = global_gp
        ## #   heatmap_legend_param = list(title_gp = title_gp, labels_gp = global_gp)
        ##     )

        #HANDLE NA values: NA goes to 0
        tmp <- as.vector(annotation[,c])
        tmp[is.na(tmp)] <- 0

        #DOUBLE ## means the param causes errors!
        foo <- HeatmapAnnotation(
            df=as.data.frame(tmp),
            ##na_col = "black",
            name = gsub("_", " ", c),
            ##col = col,            
            which='row',
            #show_legend = FALSE,
            width = unit(3, "mm"),
            ##column_names_gp = global_gp,
            ##column_title_gp = global_gp
            )
        
        ht_list <- ht_list + foo
    }
    return(ht_list)
}

heatmapSF_plot <- function(rpkmTable, annotation, plot_out, sfCorr_out) {
    #CONSTANTS
    RPKM_THRESHOLD <- 2.0
    MIN_NUM_SAMPLES_EXPRESSSING_AT_THRESHOLD <- 4
    NUM_GENES_TO_CLUSTER <- 150#250

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

    #SAVE plot
    pdf(file = plot_out, width=11,height=8.5)
    
    #make SF (sample-feature) heatmap
    Exp_data <- apply(Exp_data,1,function(x) zscore(x))
    ma_nolym <- max(Exp_data)
    mi_nolym <- min(Exp_data)
    my.breaks_nolym<-c(-3,seq(-2.5,2.5,length.out=99),3)
    param_text <- paste(RPKM_THRESHOLD, MIN_NUM_SAMPLES_EXPRESSSING_AT_THRESHOLD, NUM_GENES_TO_CLUSTER, sep=",")

    #LEN: DROP this!
    ## test<-heatmap.2(as.matrix(Exp_data),
    ##                 distfun = function(x) dist(x,method = 'euclidean'),
    ##                 hclustfun = function(x) hclust(x,method = 'ward.D2'),
    ##                 breaks=my.breaks_nolym,trace="none",scale="none",
    ##                 col=greenred(100),
    ##                 #ALLOW default of printing colum labels! labCol=F,
    ##                 cexRow=0.5,
    ##                 cexCol=0.5,
    ##                 lwid=c(0.15,0.75),
    ##                 key=FALSE, margins=c(3,7),
    ##                 #main = paste(Project_Name,"Sample-Feature Correlation"),
    ##                 main = "Sample-Feature Correlation",
    ##                 #xlab = paste(NUM_GENES_TO_CLUSTER, " genes ",Sys.Date())
    ##                 xlab = paste(NUM_GENES_TO_CLUSTER, " genes ")
    ##                 )
    
    graph2 <-Heatmap(as.matrix(Exp_data),
                     col = colorRamp2(my.breaks_nolym,  greenred(101), transparency = 0),
              #         column_dend_height = unit(2, "cm"),
                       heatmap_legend_param = list(title = "exp. level"),
                       column_title = "Sample-Feature Correlation",
                       #REMOVErow_title = "Samples",
                       show_row_names = TRUE,show_column_names = TRUE,
                       row_names_max_width = unit(3, "mm"),
                       row_names_gp = gpar(fontsize = 4),
                       column_names_gp = gpar(fontsize = 5),
                       cluster_rows = TRUE, cluster_columns=TRUE,
                     clustering_method_rows="ward.D2",
                     clustering_method_columns="ward.D2",
                     clustering_distance_rows="euclidean",
                     clustering_distance_columns="euclidean",
                     show_heatmap_legend = FALSE,
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                        )
    draw(graph2)
    
    #LEN:
    #output<-as.matrix(Exp_data)[rev(test$rowInd), test$colInd]
    output<-graph2@matrix
    #data.matrix.filename <- "XXX_bar.txt"
    #write.table(output, file=data.matrix.filename, quote=F, col.names = NA, sep="\t")
    write.table(output, file=sfCorr_out, quote=F, col.names = NA, sep="\t")
    ht_list <- Heatmap(t(as.matrix(Exp_data)),
                       col = colorRamp2(my.breaks_nolym,  greenred(101), transparency = 0),
              #         column_dend_height = unit(2, "cm"),
                       heatmap_legend_param = list(title = "exp. level"),
                       #column_title = "Sample-Feature Correlation",
                       #REMOVErow_title = "Samples",
                       show_row_names = TRUE, show_column_names = TRUE,
                       row_names_max_width = unit(3, "mm"),
                       row_names_gp = gpar(fontsize = 4),
                       #column_names_gp = gpar(fontsize = 5),
                       cluster_rows = TRUE, cluster_columns=TRUE,
                       show_heatmap_legend=FALSE,
                        )

    #LEN: not even sure if this heatmap graph is helpful!--it's super complex!
    #TURNING OFF for now!
    #ht_list <- make_complexHeatmap_annotation(ht_list, annotation)
    draw(ht_list)    

    #SAVE graphics
    dev.off()
}

args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
sf_plot_out=args[3]
sf_txt_out=args[4]

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
rownames(tmp_ann) <- tmp_ann$SampleName
samples <- intersect(colnames(rpkmTable), rownames(tmp_ann))
tmp_ann <- tmp_ann[samples,-1:-3]
#print(str(tmp_ann))
heatmapSF_plot(rpkmTable, tmp_ann, sf_plot_out, sf_txt_out)