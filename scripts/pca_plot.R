#Script to generate a PCA (principle components) plot
#Input:
#   rawRPKM_file: a comma-sep file of Gene,Sample1(RPKM),Sample2, ...,SampleN
#   annotation: annotation file in SOME TBD format
#OUTPUT:
#   plot_out: some filename ending in '.pdf' to save the plots
#   Also out plotting all as pngs in images
#BASED on garber_analysis.R by Henry Long

# load required packages
suppressMessages(library("gplots"))
suppressMessages(library("ComplexHeatmap"))
suppressMessages(library("circlize"))
suppressMessages(library("dendextend"))
suppressMessages(library("viridis"))
suppressMessages(library('dplyr'))
suppressMessages(source('viper/scripts/supp_fns.R'))

#enable stack trace
#LEN:
options(error = function() traceback(2))

pca_plot <- function(rpkmTable,annotation, RPKM_threshold,min_num_samples_expressing_at_threshold,filter_mirna,SSnumgenes, pca_plot_out) {
    
    #readin and process newdata
    newdata <- rpkmTable

    ## We want to only work with the samples that are in the meta file, so we are only selecting the count columns that are in the meta file
    newdata = newdata[,colnames(newdata) %in% rownames(tmp_ann)]    
    
    #remove genes with no RPKM values or
    newdata<-newdata[apply(newdata, 1, function(x) length(x[x>=RPKM_threshold])>min_num_samples_expressing_at_threshold),]

    #log transform of data
    newdata <- log2(newdata+1)

    ## Removing Sno and Mir mrna, parameterized
    if (filter_mirna == TRUE) {
        newdata <- newdata[ !grepl("MIR",rownames(newdata)), ]
        newdata <- newdata[ !grepl("SNO",rownames(newdata)), ]
    }

    ## Fail safe to take all genes if numgenes param is greater than what passes filters
    if (as.numeric(SSnumgenes) > nrow(newdata)) {SSnumgenes = nrow(newdata)}

    #Calculate CVs for all genes (rows)
    mean_rpkm_nolym <- apply(newdata,1,mean)
    var_rpkm_nolym <- apply(newdata,1,var)
    cv_rpkm_nolym <- abs(var_rpkm_nolym/mean_rpkm_nolym)

    #Select out the most highly variable genes into the dataframe 'Exp_data'
    Exp_data <- newdata[order(cv_rpkm_nolym,decreasing=T)[1:SSnumgenes],]

    #SAVE plot
    pdf(file = pca_plot_out)
    png_counter <- 1
    
    #Standard PCA analysis using all possible annotations
    for(c in colnames(annotation)) {
        ann <- as.matrix(annotation[, c])
        if(length(sort(unique(na.omit(as.vector(ann))))) <7) {
            ClassColors <- cmap(ann)
            
            myColors = ClassColors[ann]
            myColors[which(is.na(myColors))] <- "black"
            #png(file=paste(png_dir,"/pca_plot_",c,".png",sep=""), width = 8, height = 8, unit="in",res=300)
            png(file=paste("analysis/plots/images/pca_plot_", c, ".png", sep=""), width = 8, height = 8, unit="in",res=300) 
            png_counter <- png_counter + 1
            pca_output <- make_pca_plots(t(Exp_data), threeD = FALSE, ClassColorings = myColors, pca_title = c, legend_title =  c)
	    dev.off()
	    pca_output <- make_pca_plots(t(Exp_data), threeD = FALSE, ClassColorings = myColors, pca_title = c, legend_title =  c)
		
        }
    }

    #GET percent variances
    pc_var <- signif(100.0 * summary(pca_output)[[6]][2,], digits = 3)
    #scree plot
    #png(file=paste(png_dir,"/pca_plot_",png_counter,".png",sep=""), width = 8, height = 8, unit="in",res=300)
    png(file="analysis/plots/images/pca_plot_scree.png", width = 8, height = 8, unit="in",res=300)
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
RPKM_threshold=args[3]
min_num_samples_expressing_at_threshold=args[4]
filter_mirna = args[5]
SSnumgenes=args[6]
pca_plot_out=args[7]

#process RPKM file
# Mahesh adding check.names=F so that if there is any - or _ characters, they won't be turned to default '.'
rpkmTable <- read.table(rpkmFile, header=T, check.names=F, row.names=1, sep=",", stringsAsFactors=FALSE, dec=".")
for (n in names(rpkmTable)) {
    rpkmTable[n] <- apply(rpkmTable[n], 1, as.numeric)
}
rpkmTable = na.omit(rpkmTable)

#PROCESS ANNOTATIONS
tmp_ann <- read.delim(annotFile, sep=",", stringsAsFactors=FALSE)
#REMOVE comp_ columns
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

pca_plot(rpkmTable,tmp_ann, RPKM_threshold,min_num_samples_expressing_at_threshold,filter_mirna,SSnumgenes, pca_plot_out)
