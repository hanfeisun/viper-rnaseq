library(limma)
library(DESeq2)
library(edgeR)

limma_and_deseq_f <- function(counts, s1,s2, limma,deseq) {
    treatlist = strsplit(s2,',')[[1]]
    ctrllist = strsplit(s1,',')[[1]]

    ## Build the raw count matrix by reading in each of the count files and appending them together
    countmat <- do.call(cbind,lapply(counts,function(fn)read.table(fn,header=FALSE, sep="\t")[,3]))
    #LEN: let's take GeneName instead of GeneSymbol
    #genes <- read.table(counts[1], header=FALSE, sep="\t")[,1]
    genes <- read.table(counts[1], header=FALSE, sep="\t")[,2]
    rownames(countmat) <- genes
    
    #LEN: this doesn't work!
    #colnames(countmat) = gsub("(counts/|.read_cnt.txt)",'',counts)
    #GET SAMPLE names
    samples = sapply(counts, function(path)strsplit(path,"/")[[1]][3], simplify=FALSE, USE.NAMES=FALSE)
    #LEN: convert the above list to a vector
    samples = unlist(samples)
    colnames(countmat) = samples

    ctrllist = countmat[ ,colnames(countmat) %in% ctrllist]
    treatlist = countmat[ ,colnames(countmat) %in% treatlist]

    ntreat = ncol(treatlist)
    nctrl = ncol(ctrllist)

    data = cbind(treatlist,ctrllist)
    ##remove duplicate rownames
    #data = data[!duplicated(data),]
    
    condition = c(rep('treat',ntreat),rep('control',nctrl))
    colnames(data)=seq(ntreat+nctrl)
    preparationD <- function (countTable ,ntreat,nctrl){
        #countTable <- read.table(table, header= FALSE, row.names= 1)
        #dataPack<- data.frame(row.names= colnames(countTable), condition = c(rep('treat',ntreat),rep('control',nctrl)))
        conds <- factor(c(rep('treat',ntreat),rep('control',nctrl)))
        rownames(countTable) = make.names(rownames(countTable), unique=TRUE)
        cds <- DESeqDataSetFromMatrix(countTable, DataFrame(conds), ~ conds)
        #newCountDataSet() default:
        #newCountDataSet(countData, conditions, sizeFactors = NULL, phenoData = NULL, featureData = NULL)
        cds <- estimateSizeFactors (cds)

        #estimateSizeFactors() default:
        #estimateSizeFactors( object, locfunc=c("median", "shorth") )
        #print (sizeFactors(cds))
        #print (head(counts(cds, normalized = TRUE)))
  
        if (ntreat<=2 && nctrl<=2){
            cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
        }else{
            cds = estimateDispersions (cds)}
        
        #estimateDispersions() default:
        #estimateDispersions( object,method = c( "pooled", "pooled-CR", "per-condition", "blind" ),
        #sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
        #fitType = c("parametric", "local"),
        #locfit_extra_args=list(), lp_extra_args=list(),
        #modelFrame = NULL, modelFormula = count ~ condition, ... )
        #print(str(fitInfo(cds)))
        #################### swapped control and treat, now it gives consistent results as limma
        cds <- DESeq(cds)
        res <- results(cds)
        
        #ORDER results by XXX
        #res <- results[order(res$padj)]
        #summarize data
        #sum(res$padj < 0.1, na.rm=TRUE)

        #nbinomTest() default:
        #nbinomTest(cds, condA, condB, pvals_only = FALSE)
        return (res)
    }
    deseq_result= preparationD(data,ntreat,nctrl)

    preparationL <- function(counts,ntreat,nctrl){
        #counts <- read.delim(table, header = FALSE, row.names = 1)
        mType <- c(rep('treat',ntreat),rep('control',nctrl))
        ###In case of needing to filtration uncomment two below commands:
        #isexpr <- rowSums(cpm(counts)>1)>= each
        #counts <- counts[isexpr,]
        nf <- calcNormFactors(counts)
        #calcNormFactors() default:
        #calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"), refColumn = NULL,
        #logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
        design <- model.matrix(~mType)
        y <- voom (counts, design, lib.size = colSums(counts)*nf)
        #voom() default:
        #voom(counts, design = NULL, lib.size = NULL, normalize.method = "none", plot = FALSE, span=0.5, ...)
        fit <- lmFit(y,design)
        #lmFit() default:
        #lmFit(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation,weights=NULL,method=c("ls","robust")) 
        fit <- eBayes(fit)
        res = topTable(fit,coef=2,n=length(counts[,1]),sort="p")
        #eBayes() default:
        #eBayes(fit, proportion=0.01, stdev.coef.lim=c(0.1,4), trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
        return(res)
    }
    if (ntreat>1){
        limma_result = preparationL(data,ntreat,nctrl)
        limma_result <- cbind(ID=rownames(limma_result), limma_result)
        write.table(limma_result,limma,sep='\t',col.names=T,row.names=F,quote=F)
    }
    #LEN: setting the first column name to 'id'
    deseq_result <- cbind(id=rownames(deseq_result), as.matrix(deseq_result))
    write.table(deseq_result,deseq,sep='\t',col.names=T,row.names=F,quote=F)
}

args <- commandArgs( trailingOnly = TRUE )
arg_counts = strsplit(args[1]," ")[[1]]
#arg_s1 = strsplit(args[2], ",")
#arg_s2 = strsplit(args[3], ",")
arg_s1 = args[2]
arg_s2 = args[3]
limma_out=args[4]
deseq_out=args[5]

limma_and_deseq_f(arg_counts, arg_s1, arg_s2, limma_out, deseq_out)

        
