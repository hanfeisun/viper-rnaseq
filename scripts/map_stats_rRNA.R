library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs( trailingOnly = TRUE )

data <- read.csv( args[1], sep=",", header=TRUE )

rownames(data) <- data[,1]
data[,1] <- NULL
x <- data.frame( Sample=names(data), rRNA_Reads=as.numeric(gsub('%', '', as.matrix(data["Uniquely_mapped_reads_%",]))))
x1 <- melt(x, id.var="Sample")

png( args[2], width = 8, height = 8, unit="in",res=300 )

upper_limit <- max(x$rRNA_Reads)
limits <- seq( 0, upper_limit, length.out=10)

if( upper_limit > 1 ) { 
	limits <- round(limits)
}else {
  limits <- round(limits,digits=2)
}

colors <- c(rRNA_Reads="firebrick3")

ggplot(x1, aes(x=Sample, y=value, fill=variable)) +
geom_bar( stat = "identity" ) +
scale_y_continuous("",limits=c(0,upper_limit),  breaks=limits, labels=percent(limits/100)) +
scale_fill_manual(values=colors) +
labs( title="rRNA Report\n\n", x = "Sample Names", y="") +
guides(fill=guide_legend(title=NULL)) + theme_bw() +
#geom_hline(yintercept =5,linetype="longdash") +
#geom_hline(yintercept =10, linetype= "dash") +
theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5, size=10))

dev.off()

