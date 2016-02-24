#Script to generate a Sample-Feather heatmaps
#Input:
#   rawRPKM_file: a comma-sep file of Gene,Sample1(RPKM),Sample2, ...,SampleN
#   annotation: annotation file in SOME TBD format
#OUTPUT:
#   plot_out: some filename ending in '.pdf' to save the plots
#   sfCorr_out: a .txt file to save sample-feature correlations
#BASED on makeclustering_all_samples_HWL.R by Henry Long

# load required packages
suppressMessages(library("gplots"))
suppressMessages(library("genefilter"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("dendextend"))
suppressMessages(library("viridis"))
suppressMessages(library('dplyr'))
suppressMessages(source('snakemake/scripts/supp_fns.R'))

#enable stack trace
options(error = function() traceback(2))

heatmapSF_plot <- function(rpkmTable,tmp_ann, RPKM_threshold,min_num_samples_expressing_at_threshold,filter_mirna,SFnumgenes,num_kmeans_clust, sf_plot_out,sf_txt_out) {
    #RPKM_THRESHOLD <- 2.0
    #MIN_NUM_SAMPLES_EXPRESSSING_AT_THRESHOLD <- 4
    #NUM_GENES <- 22000 #roughly how many human genes there are
    #NUM_GENES_TO_CLUSTER <- NUM_GENES * 0.05 #cluster top 5%

    #readin and process newdata
    newdata <- rpkmTable

    #remove genes with no RPKM values or genes where not enough samples meet a minimum threshold
    newdata<-newdata[apply(newdata, 1, function(x) length(x[x>=RPKM_threshold])>min_num_samples_expressing_at_threshold),]

    #log transform of data
    newdata <- log2(newdata+1)

    ## Removing Sno and Mir mrna, parameterized
    if (filter_mirna == TRUE) {
        newdata <- newdata[-grep("MIR[[:digit:]]*",rownames(newdata)), ]
        newdata <- newdata[-grep("SNO.*",rownames(newdata)), ]
    }

    #Calculate CVs for all genes (rows)
    mean_rpkm_nolym <- apply(newdata,1,mean)
    var_rpkm_nolym <- apply(newdata,1,var)
    cv_rpkm_nolym <- abs(var_rpkm_nolym/mean_rpkm_nolym)

    #Select out the most highly variable genes into the dataframe 'Exp_data'
    Exp_data <- newdata[order(cv_rpkm_nolym,decreasing=T)[1:SFnumgenes],]

    # NOTES on clustering, not used for now
    # Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
    # Cluster options: complete (default), single, average, mcquitty, median, centroid, ward
    rowdistance = dist(as.matrix(Exp_data), method = "euclidean")
    rowcluster = hclust(rowdistance, method = "ward.D2")
    coldistance = dist(t(as.matrix(Exp_data)), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
        
    #make SF (sample-feature) heatmap
    Exp_data <- apply(Exp_data,1,function(x) zscore(x))
    ma_nolym <- max(Exp_data)
    mi_nolym <- min(Exp_data)
    my.breaks_nolym<-c(-3,seq(-2.5,2.5,length.out=99),3)
    param_text <- paste(RPKM_threshold, min_num_samples_expressing_at_threshold, SFnumgenes, sep=",")

    pdf(file = sf_plot_out, width=11,height=8.5)
    
    ha1 <- make_complexHeatmap_annotation(tmp_ann)

    if (is.numeric(num_kmeans_clust) == TRUE) {kmparam = num_kmeans_clust}
    if (is.character(num_kmeans_clust) == TRUE) {kmparam = as.numeric(unlist(strsplit(num_kmeans_clust,",")))}
    
    for (i in 1:length(kmparam)) {
        rowclusterparam = FALSE
        if (kmparam[i] == 0) {
            kmparam[i] = FALSE
            rowclusterparam = rowcluster
        }
        mapplot <- Heatmap(t(as.matrix(Exp_data)),
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = "Sample-Feature Correlation",
                     #REMOVErow_title = "Samples",
                     show_row_names = TRUE, show_column_names = TRUE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 12),
                     column_names_gp = gpar(fontsize = 8),
                     #cluster_rows = TRUE,
                     #cluster_columns=TRUE,
                     cluster_rows = rowclusterparam,
                     km = kmparam[i],
                     cluster_columns = colcluster,
                     show_heatmap_legend = FALSE,
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
        draw(mapplot)
        for(an in colnames(tmp_ann[1:ncol(tmp_ann)])) {
            decorate_annotation(an, {
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc",
                          just = "left", gp=gpar(fontsize=5), check=TRUE)
                grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc",
                          just = "right", gp=gpar(fontsize=5), check=TRUE)
            })
        }
    }
    dev.off()
        
    ## For now, repeating everything for PNG, should probably change this

    #png(file="analysis/plots/images/heatmapSF_plot.png", width = 8, height = 8, unit="in",res=300) 
    png(file="analysis/plots/images/heatmapSF_%2d_plot.png", width = 8, height = 8, unit="in",res=300)
    
    for (i in 1:length(kmparam)) {
        rowclusterparam = FALSE 
        if (kmparam[i] == 0|FALSE) {
            kmparam[i] = FALSE
            rowclusterparam = rowcluster
        }
        mapplot <-Heatmap(t(as.matrix(Exp_data)),
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = "Sample-Feature Correlation",
                     #REMOVErow_title = "Samples",
                     show_row_names = FALSE, show_column_names = TRUE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 12),
                     column_names_gp = gpar(fontsize = 8),
                     #cluster_rows = TRUE,
                     #cluster_columns=TRUE,
                     cluster_rows = rowclusterparam,
                     km = kmparam[i],
                     cluster_columns = colcluster,
                     show_heatmap_legend = FALSE,
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
        draw(mapplot)
        for(an in colnames(tmp_ann[1:ncol(tmp_ann)])) {
            decorate_annotation(an, {
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc",
                          just = "left", gp=gpar(fontsize=5), check=TRUE)
                grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc",
                          just = "right", gp=gpar(fontsize=5), check=TRUE)
            })
        }
    }
    dev.off()

    ########row_order(mapplot)
    ########column_order(mapplot)
    
    #WRITE output to file
    #output<-as.matrix(Exp_data)#[rev(Exp_data$rowInd), Exp_data$colInd]

    #if (kmparam[1] != 0) {
        #km1 = kmeans(as.matrix(t(Exp_data)), centers=kmparam[1])
        #kmclust = km1$cluster
        #output<-as.matrix(t(Exp_data))
        #out <- cbind(output, clusterNum = kmclust)

        #output<-as.matrix(t(Exp_data))
        #mapplot<-Heatmap(output, km = kmparam[1], cluster_columns = colcluster)
        #output<-output[unlist(row_order(mapplot)), unlist(column_order(mapplot))]
        #write.table(output, file=sf_txt_out, quote=F, col.names = NA, sep="\t")
    #### IN PROGRESS
    #}    

    #if (kmparam[1] == 0) {
        output<-as.matrix(t(Exp_data))
        mapplot<-Heatmap(output, cluster_rows = rowcluster, cluster_columns = colcluster)
        output<-output[unlist(row_order(mapplot)), unlist(column_order(mapplot))]
        write.table(output, file=sf_txt_out, quote=F, col.names = NA, sep="\t")
    #}
}


## Read in arguments
args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
RPKM_threshold=args[3]
min_num_samples_expressing_at_threshold=args[4]
filter_mirna=args[5]
SFnumgenes=args[6]
num_kmeans_clust=args[7]
sf_plot_out=args[8]
sf_txt_out=args[9]

## Process RPKM file
rpkmTable <- read.table(rpkmFile, check.names=F, header=T, row.names=1, sep=",", stringsAsFactors=FALSE, dec=".")
for (n in names(rpkmTable)) {
    #CONVERT to numeric!
    rpkmTable[n] <- apply(rpkmTable[n], 1, as.numeric)
}
rpkmTable = na.omit(rpkmTable)

## PROCESS ANNOTATIONS
tmp_ann <- read.delim(annotFile, sep=",", stringsAsFactors=FALSE)
## REMOVE comp_ columns
tmp_ann <- tmp_ann[ , !grepl('comp_*', names(tmp_ann))]

## Convert numerical annotations to numbers/floats
for (col in colnames(tmp_ann)) {
    ## Test first value in col for validity
    if(attr(regexpr("^\\-?\\d+\\.\\d+$",tmp_ann[1,col]), "match.length") > 0){
        #print(apply(as.matrix(tmp_ann[,col]), 2, as.numeric))
        tmp_ann[,col] <- as.vector(apply(as.matrix(tmp_ann[,col]), 2, as.numeric))
    }
}

rownames(tmp_ann) <- tmp_ann$SampleName
samples <- intersect(colnames(rpkmTable), rownames(tmp_ann))
tmp_ann <- tmp_ann[samples,-1]

## Run the function
heatmapSF_plot(rpkmTable,tmp_ann, RPKM_threshold,min_num_samples_expressing_at_threshold,filter_mirna,SFnumgenes,num_kmeans_clust, sf_plot_out,sf_txt_out)
