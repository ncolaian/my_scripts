#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(getopt)

if( length(args) != 3 ) {
  stop("Need to specify the 2 files to perform the graph on, and the outfile")
}

# 1 and 2 are the mapping stats
# Arg 3 is the outfile
#for testing
args[1] <- "/Users/ncolaian/Documents/Scott_Fungal_Work/mapping_data/res_tbl_60.txt"
args[2] <- "/Users/ncolaian/Documents/Scott_Fungal_Work/mapping_data/v2_db_mapping_stats_60.txt"

#Read in the files
stat1 <- read.delim(args[1], header=TRUE, sep="\t")
stat2 <- read.delim(args[2], header=TRUE, sep="\t")

#merge the two datasets
colnames(stat1) <- colnames(stat2)
stat1$Map_type <- "Assembly"
stat2$Map_type <- "Database"
combined_stat <- rbind(stat1,stat2)

#work when datasets are already merged
stat1$SampleID <- as.character(combined_stat$SampleID)
stat1$Map_type <- as.character(combined_stat$Map_type)
stat1$Map_type[combined_stat$Map_type == "AS"] <- "Assembly"
stat1$Map_type[combined_stat$Map_type == "DB"] <- "Genomes"
summary(stat1)

#make into data frame
stat_df <- as.data.frame(combined_stat)

stat_plot <- ggplot(stat_df, aes( x=Map_type, y=Percent) )+
  geom_boxplot(outlier.size = NA)+
  geom_jitter() +
  theme(#axis.title.x = 
        plot.title = element_text(hjust = .5),
        legend.position = "none"
        )+
  ggtitle("Assembly vs Genome Mapping (V2) at 60%")+
  labs(x="Mapping Database", y="Percentage of Total Reads Mapped")

stat_plot

ggsave("/Users/ncolaian/Documents/Scott_Fungal_Work/mapping_data/v2_db_vs_assembly_mapping.png", plot = stat_plot)
