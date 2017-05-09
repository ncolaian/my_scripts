# This Rscript will create a boxplot comparing meta analysis and CL21 bar seq

#Plot theme
barseq_theme <- function() {
  theme <- theme_replace( plot.title = element_text(hjust=.5, size = 22),
                          axis.text = element_text(size = 16),
                          axis.title = element_text(size = 16),
                          legend.title = element_text(size=18),
                          legend.text = element_text(size = 16)
  )
  theme_set(theme)
  theme_get
}

barseq_theme
library(iheatmapr)
library(plotly)
library(tidyr)
library(ggplot2)
library(plyr)

#Read in the barseq
bseq <- read.delim("/Users/ncolaian/Documents/barseq_downsample/fit_logratios_unnormalized.tab",
                   sep = "\t", header = TRUE)

#Read in the meta file
meta <- read.delim("/Users/ncolaian/Documents/Barseq_v_meta/gene_counts_id60_cog_agg_da_vec.txt",
                   sep = "\t", header = F)

#The new files will deal with the log fold changes -> use up or down but not both
lfc_meta <- read.delim("/Users/ncolaian/Documents/Barseq_v_meta/gene_counts_id60_cog_agg_tags_data.txt",
                       sep = " ")

#Format the meta and barseq data for merging
colnames(meta) <- c("Cog", "Abund")
bseq <- bseq[,c(-2,-3)]
colnames(bseq) <- c("Cog", "logr")
bseq <- unique(bseq)

#Need to merge the barseq and metafile
combined_file <- merge(bseq, meta, by = "Cog")

cor(combined_file$Abund, combined_file$logr)

#Get the top 10 COGS
bseq$Cog[order(bseq$logr, decreasing = T)[1:10]]

combined_file$Abund <- as.factor(combined_file$Abund)

#This will change the level names to what we want
combined_file$Abund <- revalue(combined_file$Abund, c("-2"="No Data", "-1"="Down", "0"="No Diff.", "1"="Up"))


# Create the box plot of log ratios and meta abundance difference
plot <- ggplot(combined_file, aes(x=Abund, y=logr))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter()+
  ggtitle(label = "Barseq LogRatio vs DAFE Abundances Day 6")+
  labs( x="DAFE Abundance Call", y="log(BK/RZ)")+
  theme( plot.title = element_text(hjust = 0.5))
plot
ggsave( "/Users/ncolaian/Documents/barseq_downsample/Barseq_v_Dafe_day6_500.png", plot)

#Perform T-tests -- This is a paired two tests 
# Null hypothesis is that there is no difference between the equal and
# the down/up
upvse <- t.test( combined_file$logr[combined_file$Abund == "Up"],
                 combined_file$logr[combined_file$Abund == "No Diff."])
dwnve <- t.test( combined_file$logr[combined_file$Abund == "Down"],
                 combined_file$logr[combined_file$Abund == "No Diff."])
upvsdown <- t.test( combined_file$logr[combined_file$Abund == "Up"],
                    combined_file$logr[combined_file$Abund == "Down"])
upvse
dwnve
upvsdown


#### LOG FOLD CHANGE WORK ####
bseq <- bseq[,c(-2,-3)]
colnames(bseq) <- c("Cog", "logr")
bseq <- unique(bseq)

#Need to merge the barseq and metafile
lf_combined_file <- merge(bseq, lfc_meta, by.x = "Cog", by.y = "row.names")

#Perform a correlation analysis
cor( lf_combined_file$logr, lf_combined_file$logFC, method = "pearson" )

#Perform the linear modeling of the data
lm_lf <- lm( logFC ~ logr, lf_combined_file )
summary(lm_lf)

lf_plot <- ggplot(lm_lf, aes(x=logFC, y=logr))+
  geom_point()+
  geom_smooth(method = lm, se=F)+
  ggtitle(label = "Barseq LogRatio vs DAFE Abundances")+
  labs( x="DAFE logFoldChange", y="Barseq LogRatio - log(BK/RZ)")+
  theme( plot.title = element_text(hjust = 0.5))
  #geom_text(data=lf_combined_file, aes(label = paste("R^2: ", r2,sep="")),parse=T,x=100,y=c(1,1.25,1.5),
            #show.legend =F)
lf_plot

ggsave( "/Users/ncolaian/Documents/barseq_downsample/90perc_scatter.png",lf_plot)

#### REMOVING OUTLIERS TRIAL ####
no_out <- lf_combined_file[lf_combined_file$logr < 6 & lf_combined_file$logr > -6,]
cor( no_out$logr, no_out$logFC, method = "pearson" )
