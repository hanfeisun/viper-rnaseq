library( ggplot2 )
library( reshape2 )

args <- commandArgs( trailingOnly = TRUE )

data <- read.table( args[1], header=TRUE )
rownames(data) <- data$Feature
sub_data <- data[ c("Introns", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons"), ]

png( args[2], width = 8, height = 8, unit="in",res=300 )

m_data <- melt( sub_data )

ggplot(data=m_data, aes(variable, value, fill=Feature)) + geom_bar(stat="identity") + xlab("Sample") + ylab("Feature") + theme_bw() + coord_flip()

