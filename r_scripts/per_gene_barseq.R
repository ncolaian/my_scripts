#This script will create a graph from the poolcount file
#that holds the number of counts for each gene


library(ggplot2)
library(tidyr) # How I will get the pool format into long format
library(reshape2)

pool <- read.delim("/Users/ncolaian/Downloads/CL21_forw_rev_genes.poolcount.txt",
                   header = TRUE, sep = "\t")
meta <- read.csv("/Users/ncolaian/Documents/CL21_barseq_metadata.csv",
                 header=TRUE)

#Wide to loong format
long_pool <- gather(pool, sample, counts, IT001:IT022)
long_pool <- merge(long_pool, meta[,c("Index.1","Sample_type")], by.x = "sample", by.y = "Index.1")

#aggregate the data by genes
agg_pool <- aggregate(pool[6:101], by=list(pool$geneid), FUN=sum)
agg_pool <- gather(agg_pool, sample, counts, IT001:IT022)
agg_pool <- merge(agg_pool, meta[,c("Index.1","Sample_type")], by.x = "sample", by.y = "Index.1")

#This graph is for the per barcode counts
graph <- ggplot(long_pool, aes(x=sample, y=counts, color = Sample_type, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Average Counts Per Barcode In Each Sample")+
  #scale_color_manual(values=c(brewer.pal(5, "Set2"), "grey60", "tan"))+
  labs( x="Sample", y= "Read Count", color="Sample Type")
ggsave(filename = "/Users/ncolaian/Documents/counts_per_bc.png", plot = graph)

#This graph is the per gene counts
gene_graph <- ggplot(agg_pool, aes(x=sample, y=log(counts+1), color = Sample_type, position_dodge(4))) +
  geom_boxplot(outlier.color = NA)+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Average Counts Per Gene in Each Barseq Sample")+
  #scale_color_manual(values=c(brewer.pal(5, "Set2"), "grey60", "tan"))+
  labs( x="Sample", y= "log(Read Count)", color="Sample Type") +
  theme( axis.text.x = element_text(size = 0),
         plot.title = element_text(hjust = 0.5))
ggsave(filename = "/Users/ncolaian/Documents/counts_gene_bseq.png", plot = gene_graph, width = 12, height = 6, units = "in")

genes <- c() 
#this for loop will get the unique gene names
for ( gene_name in pool$geneid ) {
  if ( !(gene_name %in% genes) ) {
    genes <- c(genes,gene_name)
  }
}
gene_vector <- vector(length = length(genes), mode = "integer")
names(gene_vector) <- genes
#This loop will count the number of times
for ( gene_name in pool$geneid ) {
  gene_vector[gene_name] <- gene_vector[gene_name] + 1
}
gene_vec <- gene_vector[gene_vector < 3000]

save(gene_vector, file = "/Users/ncolaian/Documents/barcode_count_vector")

