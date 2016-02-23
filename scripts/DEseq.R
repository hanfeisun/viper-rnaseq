library(limma)
library(DESeq2)
library(edgeR)

limma_and_deseq_f <- function(counts, s1,s2, limma, deseq, deseqSum_out) {
    treatlist = strsplit(s2,',')[[1]]
    ctrllist = strsplit(s1,',')[[1]]

    countmat <- read.table(counts, header=TRUE, sep=",", row.names=1)
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
        #LEN: get summary stats
        summary <- c(sum(res$padj<0.05 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.05 & res$log2FoldChange>1.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.0, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>0.5, na.rm=TRUE),
                     sum(res$padj<0.01 & res$log2FoldChange>1.0, na.rm=TRUE))
        sumTable <- matrix(summary, nrow=3, ncol=2)
        rownames(sumTable)<-c('log2FC > 0.0','log2FC > 0.5', 'log2FC > 1.0')
        colnames(sumTable)<-c('padj < 0.05','padj < 0.01')
        #LEN: write/output summary stats
        write.table(sumTable, deseqSum_out, quote=FALSE, sep=",")

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
#arg_counts = strsplit(args[1]," ")[[1]]
#arg_s1 = strsplit(args[2], ",")
#arg_s2 = strsplit(args[3], ",")
arg_counts = args[1]
arg_s1 = args[2]
arg_s2 = args[3]
limma_out=args[4]
deseq_out=args[5]
deseqSum_out=args[6]

limma_and_deseq_f(arg_counts, arg_s1, arg_s2, limma_out, deseq_out, deseqSum_out)

        
