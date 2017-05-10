#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This rcode will create bargraphs for reverse and forward files
#that has start positions on the third column

library(ggplot2)
library(getopt) # How to create a nice way to pass parameters
#library(cowplot) # How to combine 2 different plots
library(RColorBrewer)

#Passed Parameters
params= matrix(c(
  "for_count_file", "f", 1, "character",
  "rev_count_file", "r", 1, "character",
  "out_dir", "o", 1, "character"
), byrow=TRUE, ncol = 4)
opt=getopt(params)

forw_start <- read.delim(opt$for_count_file, sep=" ", header = FALSE)
rev_start <- read.delim(opt$rev_count_file, sep=" ", header = FALSE)

forw_start <- as.data.frame(forw_start)
rev_start <- as.data.frame(rev_start)

#Add column names
colnames(forw_start) <- c("Barcode", "Line", "BP")
colnames(rev_start) <- c("Barcode", "Line", "BP")

forw_graph <- ggplot(forw_start, aes(factor(BP)))+
  geom_bar()+
  ggtitle( label = "Start Positions of the Forward Anchor Sequences")+
  labs( x="Base Pair Position", y="Counts")

rev_graph <- ggplot(rev_start, aes(factor(BP)))+
  geom_bar()+
  ggtitle( label = "Start Positions of the Reverse Anchor Sequences")+
  labs( x="Base Pair Position", y="Counts")

ggsave( filename = paste(opt$out_dir, "/forward_plot.png", sep = ""), plot = forw_graph)
ggsave( filename = paste(opt$out_dir, "/reverse_plot.png", sep=""), plot = rev_graph)

