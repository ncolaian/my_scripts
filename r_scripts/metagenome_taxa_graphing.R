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
metadata <-read.delim(args[2], header = FALSE)
table.melted <- melt(table)
full_table <- merge(table.melted, metadata, by.x = "variable", by.y = "V1")
names(full_table)[names(full_table) == 'V6'] <- 'Group'
names(full_table)[names(full_table) == 'V9'] <- 'Age'
names(full_table)[names(full_table) == 'V7'] <- 'Sample_Type'
mycol <- c("variable", "Taxa_level", "Group_Name", "Group_Above", "value", "Group", "Age", "Sample_Type")
full_table <- full_table[mycol]

full_table <- as.data.frame(full_table)
levels(full_table$Group) = c(levels(full_table$Group), "BK")
full_table$Group[is.na(full_table$Group)] <- "BK"
print(summary(full_table$Group))

kingdom_subset <- full_table[full_table$Taxa_level == "K",]
phylum_subset <- full_table[full_table$Taxa_level == "P",]
class_subset <- full_table[full_table$Taxa_level == "C",]
order_subset <- full_table[full_table$Taxa_level == "O",]
family_subset <- full_table[full_table$Taxa_level == "F",]
genus_subset <- full_table[full_table$Taxa_level == "G",]
species_subset <- full_table[full_table$Taxa_level == "S",]

# Make the plots
kingdom_plot_by_sample <- ggplot(kingdom_subset, aes( x=Sample_Type, y=value)) +
  geom_boxplot(outlier.size = NA)+ #NA gets rid of multiple outliers
  geom_jitter() + #puts the points on the graph
  ggtitle("Eukaryotic Counts by Sample Type")+
  labs(x="Sample Types", y="Counts")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

kingdom_plot_byg <- ggplot(kingdom_subset, aes( x=Group, y=value, col=Group)) +
  geom_boxplot(outlier.size = NA)+ #NA gets rid of multiple outliers
  geom_point() + #puts the points on the graph
  ggtitle("Eukaryotic Counts by Genotype")+
  labs(x="Genotypes", y="Counts")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

kingdom_plot_byage <- ggplot(kingdom_subset, aes( x=Age, y=value)) +
  geom_boxplot(outlier.size = NA)+
  geom_jitter()+
  ggtitle("Eukaryotic Counts by Age")+
  labs(x="Age", y="Counts")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

#don't really need the boxplot repositioning with dodge
#Important to jitterdodge the points so they are on top of their respective boxplots
phylum_plot_byg <- ggplot(phylum_subset, aes(x=Group_Name, y=value, col = Group)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Eukaryotic Phylums by Genotype")+
  labs(x="Abundant Phylums", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

phylum_plot_byage <- ggplot(phylum_subset, aes(x=Group_Name, y=value, col = Age)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Eukaryotic Phylums by Age")+
  labs(x="Abundant Phylums", y="Counts", col="Age")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

fungi_plot_byage <- ggplot(subset(family_subset,Group_Above =="Fungi"), aes( x=Group_Name, y=value, col=Age)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Fungi Families by Age")+
  labs(x="Abundant Families", y="Counts", col="Age")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

fungi_plot_byg <- ggplot(subset(family_subset,Group_Above =="Fungi"), aes( x=Group_Name, y=value, col=Group)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Fungi Families by Genotype")+
  labs(x="Abundant Families", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

genus_plot_byg <- ggplot(subset(genus_subset,Group_Above =="Dikarya"), aes( x=Group_Name, y=value, col=Group)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Abundant Dikarya Species By Genotype")+
  labs(x="Abundant Families", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

genus_plot_byage <- ggplot(subset(genus_subset,Group_Above =="Dikarya"), aes( x=Group_Name, y=value, col=Age)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Abundant Dikarya Species By Age")+
  labs(x="Abundant Families", y="Counts", col="Age")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

genus_plot_byg_zygo <- ggplot(subset(genus_subset,Group_Above =="Zygomycota"), aes( x=Group_Name, y=value, col=Group)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Abundant Zygomycota Species By Genotype")+
  labs(x="Abundant Families", y="Counts", col="Genotype")+
  theme(plot.title = element_text(hjust = 0.5,size = 20, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

genus_plot_byage_zygo <- ggplot(subset(genus_subset,Group_Above =="Zygomycota"), aes( x=Group_Name, y=value, col=Age)) +
  geom_boxplot(outlier.size = NA, position = position_dodge(width=0.9)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.9))+
  ggtitle("Abundant Zygomycota Species By Age")+
  labs(x="Abundant Families", y="Counts", col="Age")+
  theme(plot.title = element_text(hjust = 0.5,size = 24, face="bold"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18,face = "bold"),
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16)
  )

ggsave(filename = "kingdom_plot_by_sample.png", plot = kingdom_plot_by_sample)
ggsave(filename = "kingdom_plot_byg.png", plot = kingdom_plot_byg)
ggsave(filename = "kingdom_plot_byage.png", plot = kingdom_plot_byage)
ggsave(filename = "phylum_plot_byg.png", plot = phylum_plot_byg)
ggsave(filename = "phylum_plot_byage.png", plot = phylum_plot_byage)
ggsave(filename = "fungi_plot_byg.png", plot = fungi_plot_byg)
ggsave(filename = "fungi_plot_byage.png", plot = fungi_plot_byage)
ggsave(filename = "genus_plot_byage.png", plot = genus_plot_byage)
ggsave(filename = "genus_plot_byg.png", plot = genus_plot_byg)
ggsave(filename = "zyg_genus_plot_byage.png", plot = genus_plot_byage_zygo)
ggsave(filename = "zyg_genus_plot_byg.png", plot = genus_plot_byg_zygo)
