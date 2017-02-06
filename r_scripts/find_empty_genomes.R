#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

if( length(args) != 2 ) {
  stop("Need to specify the 2 files to perform the graph on, and the outfile")
}
#For testing
#kog_table <- read.delim("/Users/ncolaian/Documents/Scott Fungal Work/metagenome_count_test_info/kog_count_tbl.txt",
 #                       header=TRUE, sep="\t", check.names = FALSE)

kog_table <- read.delim( args[1], header=TRUE, sep="\t", check.names = FALSE)


zero_table <- kog_table[,2:ncol(kog_table)][,colSums(kog_table[2:ncol(kog_table)]) == 0]

write(names(zero_table), file = args[2])
