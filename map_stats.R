library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs( trailingOnly = TRUE )

data <- read.csv( args[1], sep=",", header=TRUE )

rownames(data) <- data[,1]
data[,1] <- NULL
x <- data.frame( Sample=names(data), Total_Reads=as.numeric(as.matrix(data["Number_of_input_reads",])), Unique_Reads=as.numeric(as.matrix(data["Uniquely_mapped_reads_number",])))
x1 <- melt(x, id.var="Sample")

png( args[2], width = 8, height = 8, unit="in",res=300 )

upper_limit <- max(x$Total_Reads)
limits <- seq( 0, upper_limit, length.out=10)
colors <- c(Total_Reads="Grey", Unique_Reads="Blue")

ggplot(x1, aes(x=Sample, y=value, fill=variable)) + geom_bar( stat = "identity", position="identity") + scale_y_continuous("",limits=c(0,upper_limit), labels=comma, breaks=limits) + scale_fill_manual(values=colors) + labs( title="Read Alignment Report\n\n", x = "Sample Names", y="") + guides(fill=guide_legend(title=NULL)) + theme_bw()

dev.off()

