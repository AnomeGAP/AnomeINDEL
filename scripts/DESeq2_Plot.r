#!/usr/bin/env Rscript
library(DESeq2)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Please enter the input folder", call.=FALSE)
}

countData <- as.matrix(read.csv(paste0(args[1], "/", "gene_count_matrix.csv"), row.names="gene_id"))
database <- as.matrix(read.csv(paste0(args[1], "/", "transcript_count_matrix.csv"), row.names="transcript_id"))
condition <- factor(c("primary","primary","primary","primary","control","control","control","control","control", "control"))
#condition <- factor(strsplit(readLines(file(paste0(args[1], "/", "condition.list"), "r"), n = 1), " ")[[1]])
colData <- data.frame(row.names = colnames(database), condition)
countData <- countData[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
cat(sprintf("nrow=%d\n", nrow(dds)))
library(stringi)
dds <- dds[!(stri_detect_fixed(row.names(dds),"MSTRG."))]
cat(sprintf("nrow=%d (removing MSTRG.XXX)\n", nrow(dds)))
#dds <- dds[ rowSums(counts(dds)) > 0, ]
#cat(sprintf("nrow=%d (removing #read_count <= 1)\n", nrow(dds)))
#res <- results(dds)
#resordered <- res[order(res$padj),]
#summary(res)

rld <- rlog(dds)
colnames(rld) <- c("Parental", "Primary1", "Primary2", "Primary3", "CTC1", "CTC2", "CTC3", "Secondary1", "Secondary2", "Secondary3")
#sampleDists <- dist( t( assay(rld) ) )
sampleDists <- round( cor( assay(rld)),4 )
#write.csv(as.data.frame(resordered),file=paste0(args[1],"/", "results.csv"))
write.table(as.matrix(sampleDists),file=paste0(args[1],"/", "samples.tsv"))
write.csv(assay(rld),file=paste0(args[1],"/", "trimmed_data.csv"))
png(paste0(args[1],"/", "graph.png"))
par( mfrow = c( 2, 3))
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:3] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:3], col="#00000020", pch=20, cex=0.3 )
dev.off()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
library("ggplot2")
library("reshape2")
png(paste0(args[1],"/", "samples.png"))
sampleDistMatrix <- get_upper_tri(as.matrix( sampleDists ))
melted_sample <- melt(sampleDistMatrix, na.rm = TRUE)
#ggplot(data = melted_sample, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
ggplot(data = melted_sample, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white",
   midpoint = 0.85, limit = c(0.7,1), space = "Lab",
   name="Correlation\nCoefficient") +
  theme_minimal()+
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
dev.off()

##https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap
library(pheatmap)
library("RColorBrewer")
png(paste0(args[1],"/", "clustering.png"))
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = 'correlation', clustering_distance_cols = 'euclidean', legend_labels = 'Correlation\nCoefficient', col = colors)
#colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
#pheatmap(sampleDistMatrix, clustering_distance_rows = 'correlation', clustering_distance_cols = 'euclidean', legend_labels = 'Distance', col = colors)
dev.off()
