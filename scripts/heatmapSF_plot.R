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

    #make SF (sample-feature) heatmap
    Exp_data <- apply(Exp_data,1,function(x) zscore(x))
    ma_nolym <- max(Exp_data)
    mi_nolym <- min(Exp_data)
    my.breaks_nolym<-c(-3,seq(-2.5,2.5,length.out=99),3)
    param_text <- paste(RPKM_threshold, min_num_samples_expressing_at_threshold, SFnumgenes, sep=",")

    Exp_data = t(as.matrix(Exp_data))
    
    ha1 <- make_complexHeatmap_annotation(tmp_ann)

    coldistance = dist(t(Exp_data), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
    
    if (is.numeric(num_kmeans_clust) == TRUE) {kmparam = num_kmeans_clust}
    if (is.character(num_kmeans_clust) == TRUE) {kmparam = as.numeric(unlist(strsplit(num_kmeans_clust,",")))}
   
    pdf(file = sf_plot_out, width=11,height=8.5) 
    
    png_count = 0
    for (i in 1:length(kmparam)) {
        if (kmparam[i] == 0) {
            rowdistance = dist(Exp_data, method = "euclidean")
            rowcluster = hclust(rowdistance, method = "ward.D2")
            rowclusterparam = rowcluster
            hmdata = Exp_data
        }
        
        if (kmparam[i] != 0) {
            km1 = kmeans(Exp_data, centers=kmparam[1])
            kmclust = km1$cluster
            kmclustsort = sort(kmclust)
            ind = match(names(kmclustsort), rownames(Exp_data))
            hmdata = Exp_data[ind,]
            rowclusterparam = FALSE
        }
        
        mapplot <- Heatmap(hmdata,
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = "Sample-Feature Correlation",
                     show_row_names = TRUE, show_column_names = TRUE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 12),
                     column_names_gp = gpar(fontsize = 8),
                     cluster_rows = rowclusterparam,
                     cluster_columns = colcluster,
                     show_heatmap_legend = FALSE,
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
        
        ## First drawing into png
        png_count = png_count+1
        png(file=paste("/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/plots/images/heatmapSF_",png_count,"_plot.png",sep=""), width = 8, height = 8, unit="in",res=300)
        draw(mapplot)
        for(an in colnames(tmp_ann[1:ncol(tmp_ann)])) {
            decorate_annotation(an,
              {grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=5), check=TRUE)
              grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=5), check=TRUE)
              })
        }
        dev.off()
        
        ## Repeated to get into the pdf
        draw(mapplot)
        for(an in colnames(tmp_ann[1:ncol(tmp_ann)])) {
            decorate_annotation(an,
              {grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=5), check=TRUE)
              grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=5), check=TRUE)
              })
        }
        if (i == 1) {
            if (kmparam[1] == 0) {
                output<-Exp_data
                output<-output[unlist(row_order(mapplot)), unlist(column_order(mapplot))]
                write.table(output, file=sf_txt_out, quote=F, col.names = NA, sep="\t")
            }
            if (kmparam[1] != 0) {
                output = cbind(hmdata,kmclustsort)
                write.table(output, file=sf_txt_out, quote=F, col.names = NA, sep="\t")
            }
        }
    }
    dev.off()        
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
