#Script to generate a Sample-Sample heatmaps
#Input:
#   rawRPKM_file: a comma-sep file of Gene,Sample1(RPKM),Sample2, ...,SampleN
#   annotation: annotation file in SOME TBD format
#Output:
#   plot_out: some filename ending in '.pdf' to save the plots
#   ssSpear_out: some filename ending in 'txt' to save the sample-sample spearman correlations
#BASED on makeclustering_all_samples_HWL.R by Henry Long

## Load required packages
suppressMessages(library("gplots"))
suppressMessages(library("genefilter"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("dendextend"))
suppressMessages(library("viridis"))
suppressMessages(library('dplyr'))
suppressMessages(source('viper/scripts/supp_fns.R'))

## Enable stack trace
options(error = function() traceback(2))

heatmapSS_plot <- function(rpkmTable,tmp_ann, RPKM_threshold,min_num_samples_expressing_at_threshold,filter_mirna,SSnumgenes, ss_plot_out,ss_txt_out) {

    ## Readin and process newdata
    newdata <- rpkmTable

    ## We want to only work with the samples that are in the meta file, so we are only selecting the count columns that are in the meta file
    newdata = newdata[,colnames(newdata) %in% rownames(tmp_ann)]    
            
    ## Remove genes with no RPKM values or genes where not enough samples meet a minimum threshold
    newdata<-newdata[apply(newdata, 1, function(x) length(x[x>=RPKM_threshold])>min_num_samples_expressing_at_threshold),]

    ## Log transform of data
    newdata <- log2(newdata+1)

    ## Removing Sno and Mir mrna, parameterized
    if (filter_mirna == TRUE) {
        newdata <- newdata[ !grepl("MIR",rownames(newdata)), ]
        newdata <- newdata[ !grepl("SNO",rownames(newdata)), ]
    }

    ## Fail safe to take all genes if numgenes param is greater than what passes filters
    if (as.numeric(SSnumgenes) > nrow(newdata)) {SSnumgenes = nrow(newdata)}
    
    ## Calculate CVs for all genes (rows)
    mean_rpkm_nolym <- apply(newdata,1,mean)
    var_rpkm_nolym <- apply(newdata,1,var)
    cv_rpkm_nolym <- abs(var_rpkm_nolym/mean_rpkm_nolym)
    
    ## Select out the most highly variable genes into the dataframe 'Exp_data'
    Exp_data <- newdata[order(cv_rpkm_nolym,decreasing=T)[1:SSnumgenes],]
        
    ## Calc. spearman correlation
    cordata <- cor(Exp_data, method="spearman")

    # NOTES on clustering, not used for now
    # Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
    # Cluster options: complete (default), single, average, mcquitty, median, centroid, ward
    rowdistance = dist(as.matrix(cordata), method = "euclidean")
    rowcluster = hclust(rowdistance, method = "ward.D2")
    coldistance = dist(t(as.matrix(cordata)), method = "euclidean")
    colcluster = hclust(coldistance, method = "ward.D2")
    
    ## make SS (sample-sample) heatmap
    ma_nolym <- max(cordata)
    mi_nolym <- min(cordata)
    my.breaks_nolym<-c(mi_nolym,seq(mi_nolym + 0.01, ma_nolym - 0.01,length.out=99),ma_nolym)
    param_text <- paste(RPKM_threshold, min_num_samples_expressing_at_threshold, SSnumgenes, sep=",")

    pdf(file = ss_plot_out)

    ha1 <- make_complexHeatmap_annotation(tmp_ann)

    mapplot <-Heatmap(t(as.matrix(cordata)),
                     col = colorRamp2(my.breaks_nolym,  bluered(101), transparency = 0),
                     #column_dend_height = unit(2, "cm"),
                     #heatmap_legend_param = list(title = "exp. level"),
                     column_title = "Sample-Sample Correlation",
                     #row_title = "Samples",
                     show_row_names = TRUE, show_column_names = TRUE,
                     #row_names_max_width = unit(3, "mm"),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 6),
                     #cluster_rows = TRUE,
                     #cluster_columns=TRUE,
                     cluster_rows = rowcluster,
                     cluster_columns = colcluster,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param=list(title="corr", title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=8)),
                     #row_dend_width = unit(5, "mm"),
                     #width=unit(60,"cm"),
                     top_annotation=ha1,
                     )
    draw(mapplot)
    for(an in colnames(tmp_ann[1:ncol(tmp_ann)])) {
        decorate_annotation(an, {
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
        })
    }
    dev.off()
    
    png(file="analysis/plots/images/heatmapSS_plot.png", width = 8, height = 8, unit="in",res=300)
    draw(mapplot)
    for(an in colnames(tmp_ann[1:ncol(tmp_ann)])) {
        decorate_annotation(an, {
            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize=6), check=TRUE)
            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=6), check=TRUE)
        })
    }
    dev.off()

    #WRITE output to file
    output<-as.matrix(cordata)
    output<-output[rowcluster$order, colcluster$order]
    write.table(output, file=ss_txt_out, quote=F, col.names = NA, sep="\t")
    
}


args <- commandArgs( trailingOnly = TRUE )
rpkmFile=args[1]
annotFile=args[2]
RPKM_threshold=args[3]
min_num_samples_expressing_at_threshold=args[4]
filter_mirna = args[5]
SSnumgenes=args[6]
ss_plot_out=args[7]
ss_txt_out=args[8]

## process RPKM file
rpkmTable <- read.table(rpkmFile, header=T, row.names=1, sep=",", stringsAsFactors=FALSE, dec=".")
for (n in names(rpkmTable)) {
    # CONVERT to numeric!
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

rownames(tmp_ann) <- tmp_ann[,1]
samples <- intersect(colnames(rpkmTable), rownames(tmp_ann))
tmp_ann <- tmp_ann[samples,-1]

heatmapSS_plot(rpkmTable,tmp_ann, RPKM_threshold,min_num_samples_expressing_at_threshold,filter_mirna,SSnumgenes, ss_plot_out,ss_txt_out)
