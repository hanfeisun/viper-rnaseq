suppressMessages(library(gplots))
suppressMessages(library(VennDiagram)) 
suppressMessages(library(corrplot))

## Function used for p value correlation, used later on in the code
cor.test.p <- function(x){
     FUN <- function(x, y) cor.test(x, y)[["p.value"]]
    z <- outer(
        colnames(x),
        colnames(x),
        Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}

## Cbind function to append columsn of different lengths for the upgenes and down genes matrix
cbind.fill <- function(...){
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x)
        rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

correlation_plot_f <- function(diffiles,meta, SFnumgenes, correlation_plot,correlation_table,upvenn_plot,downvenn_plot) {

    SFnumgenes = as.numeric(SFnumgenes)

    ## NOTE: the meta file is only input to make sure you output a new correlation plot when needed

    ## Reads in a diffile list "n" long as file1,file2, etc... then do the filtering we need to and naming those file1_filt,file2_filt...
    for (i in 1:length(diffiles)) {
        tempnamestr = paste("file",i,sep="")
        #assign(tempnamestr, read.table(file=diffiles[i], col.names=c("gene","na","Gfold","EFDR","logfc","RPKM1","RPKM2"), header=FALSE))
        assign(tempnamestr, read.table(file=diffiles[i], header=TRUE, sep=","))
        
        tempnamestr2 = paste("file",i,"_filt",sep="")
        tempfile = eval(parse(text=tempnamestr))
        #assign(tempnamestr2, subset(tempfile, tempfile[,6]>1 & tempfile[,7]>1))
        assign(tempnamestr2, subset(tempfile, tempfile[,7]<0.25))
        #assign(tempnamestr2, tempfile[1:10000,])
        print(nrow(eval(parse(text=tempnamestr2))))
    }
    
    ## pulls out the gene list from each file that pass filtering for all samples
    genes_filt = intersect(file1_filt[,1], file2_filt[,1])
    if (length(diffiles) > 2) {
        for (i in 3:length(diffiles)) {
            tempnamestr = paste("file",i,"_filt",sep="")
            genes_filt = intersect(genes_filt, eval(parse(text = tempnamestr))[,1])
        }
    }
    else {
        genes_filt = intersect(file1_filt[,1], file2_filt[,1])
    }

    print(length(genes_filt))
    
    ## now with the gene list, we need to create our correlation matrix with the logfc for each gene per sample
    corrmat = data.frame(genes_filt)
    for (i in 1:length(diffiles)) {
        tempnamestr = paste("file",i,"a",sep="")
        tempfile = eval(parse(text=paste("file",i,sep="")))
        assign(tempnamestr, tempfile[tempfile[,1] %in% genes_filt,])
        print(head(file1a[,3]))
        corrmat = cbind(corrmat, eval(parse(text=paste("file",i,"a",sep="")))[,3])
    }

    #print(head(corrmat))
    
    ## gsub out the front and back for annotation purposes
    diffilenames = gsub(".deseq.csv","", diffiles)
    diffilenames = gsub("analysis/diffexp/","",diffilenames)
    colnames(corrmat)=c("gene",diffilenames)
    
    ## Found this online, it replaces the upper panels with two values, the R value or correlation, and the p value for significance. then plot it out
    panel.cor <- function(x, y, digits=2, cex.cor) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        test <- cor.test(x,y)
        Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))
        text(0.5, 0.25, cex=, paste("r=",txt))
        text(0.5, 0.75, cex=1, Signif)
    }

    ## Making sure label names fit in boxes even when you have a bunch
    diaglabel = gsub("_","\n", diffilenames)
    diaglabel = gsub("v","\nv ", diaglabel)

    ## Plot out correlation plot
    pdf(file = correlation_plot, width = 7, height = 7, pointsize = 8)
    pairs(corrmat[,2:(length(diffiles)+1)], labels = diaglabel, upper.panel=panel.cor)
    dev.off()
    
    ## Print out correlation table with samples on axes and R value in block
    rr <- abs(cor(corrmat[,1:length(diffiles)+1]))
    ptest <- cor.test.p(corrmat[,1:length(diffiles)+1])
    ptest[ptest < 0.001] <- "p<0.001"
    ptestupper = ptest[upper.tri(ptest, diag=FALSE)]
    rr[upper.tri(rr)]<-ptestupper
    rownames(rr) <- colnames(rr)
    write.csv(file = correlation_table, rr, row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    ## Part 2, create venn diagrams for comparison.

    ######## Commeneted out sections are in progress, having problem with pair function with new method
    #logfile = corrmat[,c(1,2)]
    #logfileorder = logfile[with(logfile, order(-logfile[,2])),]
    #upgenes = logfileorder[logfileorder[,2] > 1,1]
    #downgenes = logfileorder[logfileorder[,2] < -1,1]

    #if (SFnumgenes > nrow(newdata)) {SFnumgenes = nrow(newdata)}
    
    #upgenes <- data.frame(matrix(0, ncol = length(diffiles), nrow = SFnumgenes))
    upgenes <- data.frame(matrix(0, ncol = length(diffiles)))
    #upgenes <- data.frame(matrix(0, ncol = length(diffiles), nrow = length(genes_filt)))
    colnames(upgenes) = diffilenames
    #downgenes <- data.frame(matrix(0, ncol = length(diffiles), nrow = SFnumgenes))
    downgenes <- data.frame(matrix(0, ncol = length(diffiles)))
    colnames(downgenes) = diffilenames

    #print(dim(corrmat))
    
    ## 1 to 3, get rid of -1 om for statement, get rid of +1 in logfile
    for (i in 1:(length(corrmat)-1)) {
        ## Pick out the first column(genes) and the column you are sorting, order by the logfc column
        logfile = corrmat[,c(1,i+1)]
        logfileorder = logfile[with(logfile, order(-logfile[,2])),]
        #print(head(logfileorder))
        
        #print(dim(logfileorder))
        #print(head(logfileorder))
        
        #print("test")
        #upgenes[i] = logfileorder$gene[1:as.numeric(SFnumgenes)]
        #upgenes[i] = logfileorder$gene[1:SFnumgenes]
        upgenes[i] = subset(logfileorder, logfileorder[,2]>0) 
        #print("test")
        #downgenes[i] = tail(logfileorder$gene,SFnumgenes)
        downgenes[i] = subset(logfileorder, logfileorder[,2]<0)
        #downgenes[i] = tail(logfileorder$gene,as.numeric(SFnumgenes))
        #print("test")
        #upgenesnew = logfileorder[logfileorder[,2] > 1,1]
        #downgenesnew = logfileorder[logfileorder[,2] < -1,1]
        #upgenes = cbind.fill(upgenes, upgenesnew)
        #downgenes = cbind.fill(downgenes, downgenesnew)
          
    }

    #colnames(upgenes) = diffilenames
    #colnames(downgenes) = diffilenames
    #upgenes[is.na(upgenes)] <- "NOPE"
    #downgenes[is.na(downgenes)] <- "NOPE"
    #write.csv(file="/mnt/cfce-stor1/home/mgc31/code/rjlab-kallisto/gfold/output/upgenes.csv",upgenes)
    #write.csv(file="/mnt/cfce-stor1/home/mgc31/code/rjlab-kallisto/gfold/output/downgenes.csv",downgenes)

    
    panel.venn <- function(x,y) {
        par(new = TRUE)
        venn(list(x,y))
    }

    panel.overlap <- function(x,y) {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        int = intersect(x,y)
        percint = length(int)/length(x) * 100
        text(0.5,0.5,cex=1, paste("Overlap:","\n",percint,"%",sep=""))
    }

    pdf(file = upvenn_plot, width = 7, height = 7, pointsize = 8)
    pairs(upgenes, labels = diaglabel, lower.panel=panel.venn, upper.panel=panel.overlap)
    dev.off()

    pdf(downvenn_plot, width = 7, height = 7, pointsize = 8)
    pairs(downgenes, labels = diaglabel, lower.panel=panel.venn, upper.panel=panel.overlap)
    dev.off()
    
}


correlation_plot_f(
    snakemake@input[["diffiles"]],
    snakemake@input[["meta"]],
    snakemake@params[["SFnumgenes"]],
    snakemake@output[["correlation_plot"]],
    snakemake@output[["correlation_table"]],
    snakemake@output[["upvenn_plot"]],
    snakemake@output[["downvenn_plot"]]
    )

