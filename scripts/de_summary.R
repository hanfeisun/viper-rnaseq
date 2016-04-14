library(ggplot2)
library(reshape2)

args <- commandArgs( trailingOnly = TRUE )
dedata <- read.csv(args[1],header = TRUE)
df <- dedata
df[,1] <- as.character(df[,1])
mdf <- melt(df)
mdf$dir <- ifelse (grepl("up",mdf$variable),"UP","DOWN") 
mdf$adjp <- ifelse (grepl("p1",mdf$variable),"Padj < 0.1",ifelse(grepl("p05",mdf$variable),"Padj < 0.05", "Padj < 0.01"))
mdf$l2fc <- ifelse(grepl("log1",mdf$variable),"|FC| > 1", "|FC| > 2")
mdf <- subset(mdf, adjp == "Padj < 0.05" | adjp == "Padj < 0.01")
#yl <- 1.1*(max(mdf$value))
png(args[2], width = 8, height = 8, unit="in",res=300)
ggplot(mdf,aes(l2fc,value,fill=dir))+
  ggtitle("\nDifferential Gene Expression Summary\n")+
  geom_bar(position="dodge",stat="identity",width = .75)+
  scale_fill_manual("",values=c("lightblue","firebrick3"),guide = guide_legend(reverse=TRUE))+
  geom_text(aes(label=value),vjust=-0.1, position=position_dodge(width=0.9),size=3)+
  theme(text = element_text(size=10), axis.text.x = element_text(angle=0))+
  xlab("\nLog2 Fold Change\n")+
  ylab("# genes\n")+
  #coord_cartesian(ylim=c(0,yl))+
  facet_grid(adjp~Comparison,scales = "free_y")+
  theme_bw(base_size = 10)
dev.off()

