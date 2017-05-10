#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(DESeq2)
library(getopt)

params = matrix(c(
  "count_dir", "d", 1, "character",
  "genome_list", "l", 1, "character",
  "file_name", "f", 1, "character",
  "metadata", "m", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

genomes <- scan(opt$genome_list, what = character())
meta <- read.delim( opt$metadata, header = TRUE, row.names = "Id" )
meta <- as.data.frame(meta)

for( g in genomes ) {
  file_path <- paste(opt$count_dir, g, opt$file_name, sep = "/")
  count_table <- read.table(file_path, header = TRUE, sep = "\t", row.names = "name" )
  names(count_table) <- sub("X", "", names(count_table))
  count_table <- as.data.frame(count_table)
  
  
  #make sure the rows of the metadata matches the column order in the counts
  all( rownames(meta) %in% colnames(count_table) ) #make sure everything is in both
  if( !all(rownames(meta) == colnames(count_table)) ) {
    count_table <- count_table[,rownames(meta)]
  }

  dds <- DESeqDataSetFromMatrix(countData = count_table, colData = meta, design = ~ Fraction)
  
  #Think about adding this pre filtering step to get rid of stuff that doesn not have enough reads
  #dds <- dds[ rowSums(counts(dds)) > 1, ]
  
  #Specify the control
  #dds$condition <- factor(dds$condition, levels=c("BK", "RZ"))
  dds$Fraction <- relevel(dds$Fraction, ref="BK")
  
  #perform deseq and get the results
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res, row.names = rownames(res))
  
  #getting up, low, same, not enough counts - should account for every gene
  up_abund <- c( rownames(res)[res$log2FoldChange > 0 & res$padj < 0.05] )
  low_abund <- c( rownames(res)[res$log2FoldChange < 0 & res$padj < 0.05] )
  no_diff <- c( rownames(res)[res$padj >= 0.05] )
  count_issue <- c( rownames(res)[is.na(res$padj)] )
  
  print(sum(res))
  
  #Produce abundance vector
  abun_vector <- c()
  for( i in rownames(res) ) {
    if( i %in% up_abund ) {
      abun_vector <- c(abun_vector, 1)
    }
    else if ( i %in% low_abund ) {
      abun_vector <- c(abun_vector, -1)
    }
    else if ( i %in% no_diff ) {
      abun_vector <- c(abun_vector, 0)
    }
    else if ( i %in% count_issue ) {
      abun_vector <- c(abun_vector, -2)
    }
  }
  print(length(rownames(res)))
  print(length(abun_vector))
  names(abun_vector) <- rownames(res)
  
  
  #print vector
  out_f <- sub(".txt", "", file_path )
  out_f <- paste(out_f, "da_vec.txt", sep = "_")
  write.table(abun_vector, out_f, col.names = FALSE, quote = F, sep = "\t")
}



