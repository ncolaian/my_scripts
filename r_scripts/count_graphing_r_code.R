#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE) #This is needed to pass args

library(ggplot2)
library(getopt)

### Parameter variables
# The params matrix
# Each row is a parameter with 4 columns: 
#   1) long value, 
#   2) short value, 
#   3) argument type where 0 = no argument, 1 = required , 2 = optional
#   4) data type (logical, integer, double, complex, character)
# When you call this script the args are named like --verbose or -v
params = matrix(c(
  "count_file", "c", 1, "character",
  "ordered_genomes", "n", 1, "character",
  "out_file", "o", 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(params)

#tfile = "/Users/ncolaian/Documents/Scott Fungal Work/metagenome_count_test_info/kog_count_tbl_alt.txt"
#tgfile = "/Users/ncolaian/Documents/Scott Fungal Work/metagenome_count_test_info/ordered_ids.txt"

# Read in data
count_data <- read.table(opt$count_file, header = TRUE, sep="\t")
#count_data <- read.table(tfile, header = TRUE, sep="\t")

#set the order of the genomeID's
#y_label <- read.table(tgfile,header = FALSE, sep="\n")
y_label <- read.table(opt$ordered_genomes, header = FALSE, sep="\n")
y_label$V1 <- as.character(y_label$V1)
count_data$genomeID <- as.character(count_data$genomeID)

# Make data frame
count <- as.data.frame(count_data)
graphcount <- count

#make a data frame to graph where all very high points are changed to 30
graphcount$count[count$count > 30] <- 30

#make a graph
count_heat <- ggplot(graphcount, aes(x=cogID, y=genomeID)) +
  geom_tile(aes(fill=count)) +
  scale_fill_gradientn(colours = colorRampPalette(c("yellow", "red", "darkred"))(20),
                       values = c(0, .3, 1), name = "Copy Number") +
  scale_y_discrete(limits = y_label$V1) +
  theme( axis.text.y = element_blank(), 
         axis.ticks.y = element_blank(),
         #axis.text.x = element_text( angle=90 ),
         axis.text.x = element_blank(), #for full dataset
         plot.title = element_text(hjust = 0.5) ) +
  ggtitle( paste("COG Copy Number Per Genome", sep="") )+
  labs( x="COG", y = "Genomes" )

#count_heat
#save the graph to the out_file
ggsave(opt$out_file, plot = count_heat, device = NULL)
