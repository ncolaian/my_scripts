#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #This is needed to pass args

library(ggplot2)
#Pass in the count table file and then the metadata file
if ( length(args)!=3 ) {
  stop("Need to specify the csv taxa count file, and the metadata file with genotype and age data")
}

count_table <- read.delim(args[1], header=TRUE)
metadata <- read.delim(args[2], header=FALSE)

#This line merges the table based on the id's.
full_table <- merge(count_table,metadata, by.x="Id", by.y = "V1")

#Changing the names on the metadata table
names(full_table)[names(full_table) == "V6"] <- "Genotype"
names(full_table)[names(full_table) == "V7"] <- "Fraction"
names(full_table)[names(full_table) == "V9"] <- "Age"

#Get a subset of the full table 
#1) Create a vector 
#2) Pull out the parts of the matrix with the vector
col_wanted <- c( "Id", "Perc_map", "Genotype", "Fraction", "Age" )
full_table <- full_table[col_wanted]

#Make the table a dataframe
full_table <- as.data.frame(full_table)

#Need to make NA show up ( Kinda forgot how or why this happened )
levels(full_table$Genotype) = c(levels(full_table$Genotype), "NA")
full_table$Genotype[is.na(full_table$Genotype)] <- "NA"
print(summary(full_table$Genotype)) # Check to make sure we are getting correct results

#Produce the graph using ggplot2
mapped_back_graph <- ggplot( full_table, aes( x=Fraction, y=Perc_map, col=Genotype))+
  geom_boxplot(outlier.size = NA,position = position_dodge(width=0.9),fill="white")+
  geom_point(position=position_jitterdodge(dodge.width = 0.9), aes(col=Genotype, shape=Age), size = 2.5)+
  ggtitle(paste(args[3], "% Identity Mapping to Scaffolds", sep="" ))+
  labs( x= "Fraction", y= "Percent of Reads Mapped")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

#save plot 
ggsave( filename = paste("Percent_mapped_graph_diff_",args[3],".png", sep="" ), plot = mapped_back_graph)



