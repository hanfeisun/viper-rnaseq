args <- commandArgs( trailingOnly = TRUE )

data <- read.csv( args[1], sep=",", header=TRUE )

rownames(data) <- data[,1]

data[,1] <- NULL

data_matrix <- as.matrix( data[c("Number_of_input_reads", "Uniquely_mapped_reads_number"),] )

png( args[2] )

#par(mar=c(max(5.1,max(nchar(colnames(data_matrix)))/1.8),12.1,4.1 ,2.1))
par(mar=c(5.1,15.1,4.1,2.1))  

barplot( as.numeric(data_matrix["Number_of_input_reads",]), border=TRUE, main="Mapping Summary", las=2, names.arg=colnames(data), axes=FALSE, ylim=c(0,max(as.numeric( data_matrix["Number_of_input_reads",] ))) )

par(new=TRUE)

barplot( as.numeric(data_matrix["Uniquely_mapped_reads_number",]), col="blue", axes=TRUE, ylim=c(0,max(as.numeric( data_matrix["Number_of_input_reads",] ) ) ) )


par(xpd=TRUE, oma=c(1,5,1,1))
legend( -(length(colnames(data_matrix))/0.80),max(as.numeric( data_matrix["Number_of_input_reads",] ) )/2, c("Total Raw", "Unique Mapped"), cex=0.9, fill=c("grey", "blue"))

