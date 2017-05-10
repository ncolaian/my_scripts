# This script will look at the log ratio -> from paper = fitness

#Function that basically bins the values


log_ratio <- read.delim("/Users/ncolaian/Documents/barseq_bygroup/fit_logratios_unnormalized.tab",
           sep = "\t", header = TRUE)
meta <- read.delim("/Users/ncolaian/Documents/FEBA_comb_samp_metadata.txt",
                   header=T, sep = "\t")
log_ratio$desc <- NULL
log_ratio <- unique(log_ratio)

new_headers <- as.character(unique(meta$Index))[-3]
colnames(log_ratio)[3:12] <- new_headers

# The 6,12,15 samples are 7,9,11 ( BS, RZ, EC )
BS_vs_RZ <- log_ratio$IT009 - log_ratio$IT007 #RZ - BS
BS_vs_EC <- log_ratio$IT011 - log_ratio$IT007 #EC - BS
RZ_vs_EC <- log_ratio$IT011 - log_ratio$IT009 #EC - RZ

# Add row names
names(BS_vs_RZ) <- log_ratio$locusId
names(BS_vs_EC) <- log_ratio$locusId
names(RZ_vs_EC) <- log_ratio$locusId

BS_vs_RZ <- round(BS_vs_RZ)
#Load ggplot

histogram(BS_vs_RZ, type="count", xlab = "Rounded(Log_ratio(RZ) - Log_ratio(BS))")
histogram(BS_vs_EC, type="count", xlab = "Rounded(Log_ratio(EC) - Log_ratio(BS))")
histogram(RZ_vs_EC, type="count", xlab = "Rounded(Log_ratio(EC) - Log_ratio(RZ))")

#get the clustering done
hclust_log_ratio <- log_ratio[,-(1:2)]
hclust_log_ratio <- sapply(hclust_log_ratio, as.numeric)
rownames(hclust_log_ratio) <- log_ratio$locusId
col_hclust <- hclust(as.dist(1-cor(hclust_log_ratio)), "ward.D")

# Get metadata ready for table
meta <- meta[meta$Index != "IT003",]
condense_meta <- meta[,c("Index", "Day_harvest", "Sample_type")]
condense_meta$Day_harvest <- as.character(condense_meta$Day_harvest)
condense_meta$Sample_type <- as.character(condense_meta$Sample_type)
condense_meta[is.na(condense_meta)] <- "N/A"


#Heatmap creation
library(iheatmapr)
library(plotly)

my_heatmap <- main_heatmap(hclust_log_ratio, colors=c("Blue", "Light Yellow", "Red"), name="Log_Ratios vs. Inoculum" ) %>%
  add_col_annotation(condense_meta[c("Day_harvest", "Sample_type")]) %>%
  add_col_dendro(col_hclust) %>%
  add_col_title( "Combined Samples (n=10)", side="bottom", font=list(size=16)) %>%
  add_row_title("Cogs", side="left", font=list(size=16)) %>%
  add_col_title( "Comparisons of the Log Ratios of Mean Inoculum vs Comb. Samples", side="top", font=list(size=16))
my_heatmap

my_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/barseq_bygroup/log_ratio_comparison.png")


#### This is for the long day stuff ####
#get the clustering done
hclust_log_ratio <- log_ratio[,c(8,10,12)]
hclust_log_ratio <- sapply(hclust_log_ratio, as.numeric)
rownames(hclust_log_ratio) <- log_ratio$locusId
col_hclust <- hclust(as.dist(1-cor(hclust_log_ratio)), "ward.D")

# Get metadata ready for table
meta <- meta[meta$Index != "IT003",]
condense_meta <- meta[,c("Index", "Day_harvest", "Sample_type")]
condense_meta$Day_harvest <- as.character(condense_meta$Day_harvest)
condense_meta$Sample_type <- as.character(condense_meta$Sample_type)
late_day_meta <- condense_meta[c(6,8,10),]

my_heatmap <- main_heatmap(hclust_log_ratio, colors=c("Blue", "Light Yellow", "Red"), name="Log_Ratios vs. Inoculum" ) %>%
  add_col_annotation(late_day_meta["Sample_type"]) %>%
  add_col_dendro(col_hclust) %>%
  add_col_title( "Combined Samples (n=10)", side="bottom", font=list(size=16)) %>%
  add_row_title("Cogs", side="left", font=list(size=16)) %>%
  add_col_title( "Long Day Log Ratios of Mean Inoculum vs Comb. Samples", side="top", font=list(size=16))
my_heatmap

my_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/barseq_bygroup/ld_log_ratio_comparison.png")
