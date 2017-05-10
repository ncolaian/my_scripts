library(ggplot2)
library(cowplot) # How to combine 2 different plots
library(RColorBrewer)
library(tidyr)

barseq_theme <- function() {
  theme <- theme_replace( plot.title = element_text(hjust=.5, size = 22),
                 axis.text = element_text(size = 16),
                 axis.title = element_text(size = 16),
                 legend.title = element_text(size=16),
                 legend.text = element_text(size = 14)
  )
  theme_set(theme)
  theme_get
}
barseq_theme()
#####
#This r code graphs the metadata for the barseq experiments.
#It graphs the reads used and total reads
#There is also an example on how to combine graphs together into a single figure
#This file also creates a graph comparing 1 mismatch

orig_meta <- read.csv("/Users/ncolaian/Documents/CL21_barseq_metadata.csv",
                      header=TRUE)

count_data <- read.csv("/Users/ncolaian/Downloads/barseq_usage_values.csv",
                       header=TRUE)

new_pool <- read.delim("/Users/ncolaian/Downloads/CL21_for_rev_only_genes_poolcount.txt",
                       sep = "\t", header = TRUE)

orig_meta <- as.data.frame(orig_meta)
#count_data <- as.data.frame(new_pool)

######SKIP TO THE PORTION CODED FOR NEW DATA########

#combined <- merge.data.frame(orig_meta, count_data, by.x = "Index.1", by.y = "Index")

#I am changing the levels of the factor so that the bargraphs will be graphed in order
levels(orig_meta$Day_harvest) <- c("0", "6", "12", "15", "Con", "UK", "NA")

#This is the plot to look at the different treatments For Usable Barcodes
used_reads <- ggplot(orig_meta, aes(x=Sample_type, y=Read.Counts, color = Day_harvest, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Total Usable Barseq Reads for CL21")+
  labs( x="Sample Type", y= "Read Count", color="Day Harvested")+
  scale_color_manual(values=c(brewer.pal(5, "Set2"), "grey60", "tan"))+
  theme( plot.title = element_text(hjust = .5), 
         axis.text.x = element_text(size = 9),
         legend.title = element_text( size = 11),
         legend.text = element_text( size = 9) )+
  geom_hline(yintercept=1500000, colour = "grey30", linetype="dashed")

all_reads <- ggplot(orig_meta, aes(x=Sample_type, y=Read.Counts, color = Day_harvest, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Total Reads for CL21 Barseq")+
  labs( x="Sample Type", y= "Read Count", color="Day Harvested")+
  scale_color_manual(values=c(brewer.pal(5, "Set2"), "grey60", "tan"))+
  theme( plot.title = element_text(hjust = .5), 
         axis.text.x = element_text(size = 9),
         legend.title = element_text( size = 11),
         legend.text = element_text( size = 9) )+
  geom_hline(yintercept=2000000, colour = "grey30", linetype="dashed")

used_reads_day <- ggplot(orig_meta, aes(x=Day_harvest, y=nUsed, color = Sample_type, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Total Usable Barcode Reads for CL21")+
  scale_color_manual(values=c(brewer.pal(8, "Set2"), "darkred"))+
  labs( x="Day Harvested", y= "Read Count", color="Sample Type")+
  theme( plot.title = element_text(hjust = .5), 
         axis.text.x = element_text(size = 9),
         legend.title = element_text( size = 11),
         legend.text = element_text( size = 9) )+
  geom_hline(yintercept=1500000, colour = "grey30", linetype="dashed")


all_reads
used_reads
used_reads_day

#This is how you can combine 2 different graphs onto the same file
combined$full_read_dif <- combined$Read.Counts - combined$nReads
combined$used_read_dif <- combined$nReads - combined$nUsed

diff_plot <- ggplot(combined, aes(x=Well, y=full_read_dif)) +
  geom_boxplot(outlier.color = NA)+
  geom_jitter()+
  scale_y_continuous(limits = c(0, 7000000))+
  ggtitle(label = "Read.count - nReads")+
  labs( x="Well", y= "Total Reads Difference", color="Rep")+
  theme( plot.title = element_text(hjust = .5) )

used_diff <- ggplot(combined, aes(x=Well, y=used_read_dif)) +
  geom_boxplot(outlier.color = NA)+
  geom_jitter()+
  scale_y_continuous(limits = c(0, 7000000))+
  ggtitle(label = "nReads - nUsed")+
  labs( x="Well", y= "Total Reads Difference")+
  theme( plot.title = element_text(hjust = .5) )

plot_grid(diff_plot, used_diff,
          ncol = 2, nrow = 1)

df = data.frame(diff = c(combined$full_read_dif, combined$used_read_dif))
df$group = c(rep("Primer", length(combined$full_read_dif)), rep("Barcode", length(combined$used_read_dif)))


ggplot(df, aes(x=group, y=diff)) + 
  geom_boxplot(outlier.color = NA) +
  geom_jitter()+
  ggtitle( label = "Amount of Reads Not Used During Barseq") +
  labs ( x = "Why Reads Were Not Used", y="Reads Dropped") +
  theme( plot.title = element_text(hjust=.5))

### Graphing the different 
mismatched <- read.delim("/Users/ncolaian/Downloads/CL21_forward_2.txt", header = TRUE,
           sep = "\t")
regular <- read.delim("/Users/ncolaian/Downloads/CL21_forward.txt", header=TRUE,
                      sep="\t")
mismatched$nUsed <- mismatched$nUsed - regular$nUsed
mismatched$nReads <- regular$nReads - regular$nReads

ggplot(mismatched, aes(x=Index,y=nReads), stat="identity")+
  geom_boxplot()+
  ggtitle( label = "Read Difference Between Mismatch or Not")+
  theme( plot.title = element_text(hjust=.5, size = 11))


#Create a boxplot with the three different ways reads drop out
####SHOULD START HERE FOR NEW DATA####

col_sum <- colSums(new_pool[,6:101])
col_sum <- as.data.frame(col_sum)
col_sum$Index <- colnames(new_pool)[6:101]
count_data <- merge(col_sum, count_data, by = "Index" )
count_data <- merge(count_data, orig_meta[,c("Index.1", "Read.Counts", "Sample_type", "Day_harvest")], by.x = "Index", by.y = "Index.1")
lossed_reads_df <- data.frame(count_data$Index, (count_data$Read.Counts-count_data$nReads),
                              (count_data$nReads-count_data$nUsed), (count_data$nUsed-count_data$col_sum))
colnames(lossed_reads_df) <- c("Index", "No Anchor", "Bad Barcode", "Intergenic Barcode")
lossed_reads_long <- gather(lossed_reads_df, LossType, Counts, -Index)
lossed_reads_long$Counts[lossed_reads_long$Counts < 0] <- 0 
lossed_reads_long$LossType <- factor(lossed_reads_long$LossType, levels = c("No Anchor", "Bad Barcode", "Intergenic Barcode"))
lossed_reads_long <- merge(lossed_reads_long, count_data[,c("Index", "Sample_type", "Day_harvest")], by="Index")
lossed_reads_long$Day_harvest <- ordered(lossed_reads_long$Day_harvest, c("0", "6", "12", "15", "Con", "UK", "N/A"))
lossed_reads_long$Sample_type <- ordered(lossed_reads_long$Sample_type,
                                         c("B_LIBINOC", "A_LIBINOC", "BLANK", "N/A", "BULK_SOIL", "NEIGH_SOIL", "RHIZ", "RINSE_ROOT", "EC"))

blues <- brewer.pal(n=5, "Blues")

#for day_harvest
plot <- ggplot(lossed_reads_long, aes(x=LossType, y=Counts))+
  geom_boxplot(outlier.size = NA)+
  geom_jitter(aes(color=lossed_reads_long$Day_harvest), position = position_jitterdodge())+
  ggtitle( label = "Number of Reads Excluded in CL21 Barseq")+
  labs( x = "Type of Loss", y = "Total Reads Lost", color="Day Harvest")+
  scale_color_manual(values=c(blues[2], blues[3], blues[4], blues[5], "black", "tan", "purple"))+
  theme( plot.title = element_text(hjust = .5))
plot

#for sample type
plot_samp <- ggplot(lossed_reads_long, aes(x=LossType, y=Counts))+
  geom_boxplot(outlier.size = NA)+
  geom_jitter(aes(color=lossed_reads_long$Sample_type), position = position_jitterdodge())+
  ggtitle( label = "Number of Reads Excluded in CL21 Barseq")+
  labs( x = "Type of Loss", y = "Total Reads Lost", color="Sample_Type")+
  scale_color_manual(values=c("darkgray", "lightgray", "black","lightblue", "saddlebrown", "sandybrown","gold", "lightgreen", "green"))+
  theme( plot.title = element_text(hjust = .5))

ggsave("/Users/ncolaian/Documents/barseq_reads_dropped_day.png", plot, width = 9)
ggsave("/Users/ncolaian/Documents/barseq_reads_dropped_samp.png", plot_samp, width = 9)

orig_meta <- merge(orig_meta, col_sum, by.x = "Index.1", by.y = "Index")

orig_meta$Day_harvest <- ordered(orig_meta$Day_harvest, c("0", "6", "12", "15", "Con", "UK", "N/A"))

#This is the plot to look at the different treatments For Usable Barcodes
used_reads <- ggplot(orig_meta, aes(x=Sample_type, y=col_sum, color = Day_harvest, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Total Usable Barseq Reads for CL21")+
  labs( x="Sample Type", y= "Read Count", color="Day Harvested")+
  scale_color_manual(values=c(blues[2], blues[3], blues[4], blues[5], "black", "tan", "purple"))+
  theme( plot.title = element_text(hjust = .5))+
  geom_hline(yintercept=1500000, colour = "grey30", linetype="dashed")

all_reads <- ggplot(orig_meta, aes(x=Sample_type, y=Read.Counts, color = Day_harvest, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Total Reads for CL21 Barseq")+
  labs( x="Sample Type", y= "Read Count", color="Day Harvested")+
  scale_color_manual(values=c(blues[2], blues[3], blues[4], blues[5], "black", "tan", "purple"))+
  theme( plot.title = element_text(hjust = .5))+
  geom_hline(yintercept=2000000, colour = "grey30", linetype="dashed")

used_reads_day <- ggplot(orig_meta, aes(x=Day_harvest, y=col_sum, color = Sample_type, position_dodge(4))) +
  geom_boxplot(outlier.color = NA, position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge())+
  ggtitle(label = "Total Usable Barcode Reads for CL21")+
  scale_color_manual(values=c(brewer.pal(8, "Set2"), "darkred"))+
  labs( x="Day Harvested", y= "Read Count", color="Sample Type")+
  theme( plot.title = element_text(hjust = .5))+
  geom_hline(yintercept=1500000, colour = "grey30", linetype="dashed")

barseq_theme()
ggsave("/Users/ncolaian/Documents/total_reads.png", all_reads, width = 11)
ggsave("/Users/ncolaian/Documents/total_usable_reads.png", used_reads, width = 11)
ggsave("/Users/ncolaian/Documents/day_usable_reads.png", used_reads_day, width = 11)
used_reads
used_reads_day


