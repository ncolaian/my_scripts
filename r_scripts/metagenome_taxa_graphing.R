#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) #This is needed to pass args

library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)
#Pass in the table file and then the metadata file
if ( length(args)!=2 ) {
  stop("Need to specify the csv taxa count file, and the metadata file with genotype and age data")
}
  
table <- read.csv(args[1], header = TRUE, check.names=F)
metadata <-read.delim(args[2], header = FALSE )
table.melted <- melt(table)
full_table <- merge(table.melted, metadata, by.x = "variable", by.y = "V1")
names(full_table)[names(full_table) == 'V6'] <- 'Group'
names(full_table)[names(full_table) == 'V9'] <- 'Age'
mycol <- c("variable", "Taxa_level", "Group_Name", "Group_Above", "value", "Group", "Age")
full_table <- full_table[mycol]

full_table <- as.data.frame(full_table)

kingdom_subset <- full_table[full_table$Taxa_level == "K",]
phylum_subset <- full_table[full_table$Taxa_level == "P",]
class_subset <- full_table[full_table$Taxa_level == "C",]
order_subset <- full_table[full_table$Taxa_level == "O",]
family_subset <- full_table[full_table$Taxa_level == "F",]
genus_subset <- full_table[full_table$Taxa_level == "G",]
species_subset <- full_table[full_table$Taxa_level == "S",]

# Make the plots
kingdom_plot_byg <- ggplot(kingdom_subset, aes( x=Group, y=value, col=Group)) +
  geom_boxplot(outlier.size = NA)+ #NA gets rid of multiple outliers
  geom_point() + #puts the points on the graph
  ggtitle("Eukaryotic Counts by Genotype")+
  labs(x="Genotypes", y="Counts")+
  theme(plot.title = element_text(hjust = 0.5))

kingdom_plot_byage <- ggplot(kingdom_subset, aes( x=Age, y=value)) +
  geom_boxplot(outlier.size = NA)+
  geom_point()+
  ggtitle("Eukaryotic Counts by Age")+
  labs(x="Age", y="Counts")+
  theme(plot.title = element_text(hjust = 0.5))

#don't really need the boxplot repositioning with dodge
#Important to jitterdodge the points so they are on top of their respective boxplots
phylum_plot_byg <- ggplot(phylum_subset, aes(x=Group_Name, y=value, col = Group)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Eukaryotic Phylums by Genotype")+
  labs(x="Significant Phylums", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5))

phylum_plot_byage <- ggplot(phylum_subset, aes(x=Group_Name, y=value, col = Age)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Eukaryotic Phylums by Age")+
  labs(x="Age", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5))

fungi_plot_byage <- ggplot(subset(family_subset,Group_Above =="Fungi"), aes( x=Group_Name, y=value, col=Age)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Fungi Families by Age")+
  labs(x="Age", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5))

fungi_plot_byg <- ggplot(subset(family_subset,Group_Above =="Fungi"), aes( x=Group_Name, y=value, col=Group)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Fungi Families by Genotype")+
  labs(x="Significant Families", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = "kingdom_plot_byg.png", plot = kingdom_plot_byg)
ggsave(filename = "kingdom_plot_byage.png", plot = kingdom_plot_byage)
ggsave(filename = "phylum_plot_byg.png", plot = phylum_plot_byg)
ggsave(filename = "phylum_plot_byage.png", plot = phylum_plot_byage)
ggsave(filename = "fungi_plot_byg.png", plot = fungi_plot_byg)
ggsave(filename = "fungi_plot_byage.png", plot = fungi_plot_byage)