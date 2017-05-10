### My first attempt at a complex heatmap

#This is the theme we will use to creat ggplot graphs now
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

randomize_func <- function(x) {
  num_of_func <- nchar(x)
  pos <- sample(1:num_of_func, 1)
  new_func<- substr(x, pos, pos)
  return(new_func)
}

transform_matrix_to_tile <- function(x) {
  x <- as.matrix(x)
  x[x<100] <- 0
  x[(x>=100 & x<200)] <- 1
  x[(x>=200 & x<300)] <- 2
  x[x>=300] <- 3
  return(x)
}

library(ggplot2)
library(tidyr)
library(reshape2)
library(iheatmapr)
library(plotly)
library(datasets)
library(vegan)

barseq_theme()

annote <- read.delim("/Users/ncolaian/Documents/barseq_bygroup/all_annote.txt",
                     header = T, sep="\t")
cognames <- read.delim("/Users/ncolaian/Documents/barseq_bygroup/cognames2003-2014.tab",
                       header=T, sep="\t")
pool <- read.delim("/Users/ncolaian/Documents/barseq_bygroup/CL21_for_genes.poolcount",
                   header=T, sep="\t")
meta <- read.delim("/Users/ncolaian/Documents/FEBA_comb_samp_metadata.txt",
                   header=T, sep = "\t")
counts_used <- c(10811921, 10689888, 159, 8775461, 10157993, 51392015, 39058982, 4239612, 9607230, 22829732, 10637071)


#Need to put the put a different record for each IT
long_pool <- gather(pool, sample, counts, IT001:IT011)

#need to merge all the data
long_pool <- merge(long_pool, meta[,c("Index","Sample_type", "Day_harvest")], by.x = "sample", by.y = "Index")
long_pool_cog <- merge(long_pool, annote[,c("cog", "locus_tag")], by.x = "geneid", by.y = "locus_tag")
long_pool_func <- merge(long_pool_cog, cognames[,c("COG", "func")], by.x = "cog", by.y = "COG")

#For the funcs you need to randomly assign it a one var. func
long_pool_func$func <- sapply(as.character(long_pool_func$func), randomize_func)

#aggregate based on Cluster
lust_long_pool_cog <- aggregate(long_pool_cog$counts,by=list(long_pool_cog$cog, long_pool_cog$sample), FUN="sum")
colnames(lust_long_pool_cog) <- c("Cog","Samples", "Counts")
lust_long_pool_cog <- spread(lust_long_pool_cog, Samples, Counts)
row.names(lust_long_pool_cog) <- as.character(lust_long_pool_cog$Cog)
lust_long_pool_cog$Cog <- NULL

#aggregate based on function
lust_long_pool_func <- aggregate(long_pool_func$counts, by=list(
  long_pool_func$func, long_pool_func$sample), FUN="sum")
colnames(lust_long_pool_func) <- c("Function", "Samples", "Counts")
lust_long_pool_func <- spread(lust_long_pool_func, Samples, Counts)
row.names(lust_long_pool_func) <- lust_long_pool_func$Function
lust_long_pool_func$Function <- NULL

#Normalize the counts -> old normalization
#names(counts_used) <- colnames(lust_long_pool_cog)
#counts_used <- counts_used/mean(counts_used)
#take out bad sample 
lust_long_pool_cog$IT003 <- NULL
lust_long_pool_func$IT003 <- NULL

#Rarefy to the minimum count. The rrarefy works on rows so you have to transform
#the matrix and then transform it back
min = sort(apply(lust_long_pool_cog, 2, sum))[1]
lust_long_pool_cog <- t(rrarefy(t(lust_long_pool_cog), min))
lust_long_pool_func <- t(rrarefy(t(lust_long_pool_func), min))
tile_plot_cog <- lust_long_pool_cog
tile_plot_func <- lust_long_pool_func


#heatmap prep
#1) the counts need to be log10
lust_long_pool_cog <- log10(as.matrix(lust_long_pool_cog+1))
lust_long_pool_func <- log10(as.matrix(lust_long_pool_func+1))

#tilemap prep
tile_plot_cog <- transform_matrix_to_tile(tile_plot_cog)
tile_plot_func <- transform_matrix_to_tile(tile_plot_func)

#2)Cluster analysis --> Creates a distance matrix and then computes thedistance btwn samples
col_hclust <- hclust(as.dist(1-cor(lust_long_pool_cog)), "ward.D")
col_hclust_func <- hclust(as.dist(1-cor(lust_long_pool_func)), "ward.D")
tile_hclust_cog <- hclust(as.dist(1-cor(tile_plot_cog)), "ward.D")
tile_gclust_func <- hclust(as.dist(1-cor(tile_plot_func)), "ward.D")

##### MAKE COG HEATMAP ########

### SETUP ###
#get the sample types
type = meta[,c("Index", "Sample_type")]

#Get the time values
time = meta[,c("Index", "Day_harvest")]
time$Day_harvest <- as.character(time$Day_harvest)
time[is.na(time)] <- "N/A"
type$Time <- time$Day_harvest

type <- type[-3,] # this got rid of the bad inoculum sample

### CREATE HEATMAP FOR COGS

my_heatmap <- main_heatmap(lust_long_pool_cog, colors="YlOrRd", 
                           name="Log10(counts)") %>%
  add_col_annotation(type[,c("Sample_type", "Time")]) %>%
  add_col_dendro(col_hclust) %>%
  add_col_title("Combined Samples ( n=10 )", side="bottom", font=list(size=16)) %>%
  add_row_title("Cogs", side="left", font=list(size=16)) %>%
  add_col_title("Rarefied Cog Counts from Combined Barseq Samples", side="top", font=list(size=18))
my_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/cog_comb_counts_rare.png")
my_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/cog_comb_counts_rare.html")

my_heatmap
#### Make the Function Heatmap ####
func_heatmap <- main_heatmap(lust_long_pool_func, colors="YlOrRd",
                             name = "Log10(counts)") %>%
  add_row_labels() %>%
  add_col_annotation(type[,c("Sample_type", "Time")]) %>%
  add_col_dendro(col_hclust_func) %>%
  add_col_title("Combined Samples ( n=10 )", side="bottom", font=list(size=16)) %>%
  add_row_title("Functions", side="left", font=list(size=16)) %>%
  add_col_title("Rarefied Cog Fxn Counts from Combined Barseq Samples", side="top", font=list(size=18))


func_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/fxn_comb_counts_rare.png")
func_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/fxn_comb_counts_rare.html")

func_heatmap

##### MAKE TILE PLOT #####

#Cog 1st
cog_tile_hm <- main_heatmap(tile_plot_cog, colors=c("white","red","blue","green"),
                             name = "Tile Scale") %>%
  #add_row_labels() %>%
  add_col_annotation(type[,c("Sample_type", "Time")]) %>%
  add_col_dendro(tile_hclust_cog) %>%
  add_col_title("Combined Samples ( n=10 )", side="bottom", font=list(size=16)) %>%
  add_row_title("Cogs", side="left", font=list(size=16)) %>%
  add_col_title("Tiled Cog Counts from Combined Barseq Samples", side="top", font=list(size=18))

### Function plot
func_tile_hm <- main_heatmap(tile_plot_func, colors=c("white","red","blue","green"),
                            name = "Tile Scale") %>%
  add_row_labels() %>%
  add_col_annotation(type[,c("Sample_type", "Time")]) %>%
  add_col_dendro(tile_hclust_cog) %>%
  add_col_title("Combined Samples ( n=10 )", side="bottom", font=list(size=16)) %>%
  add_row_title("Functions", side="left", font=list(size=16)) %>%
  add_col_title("Tiled Cog Fxn Counts from Combined Barseq Samples", side="top", font=list(size=18))

cog_tile_hm
func_tile_hm

func_tile_hm %>% save_iheatmap("/Users/ncolaian/Documents/barseq_bygroup/tiled_rare_counts_func.png")
func_tile_hm %>% save_iheatmap("/Users/ncolaian/Documents/barseq_bygroup/tiled_rare_counts_func.html")
cog_tile_hm %>% save_iheatmap("/Users/ncolaian/Documents/barseq_bygroup/tiled_rare_counts_cog.png")
cog_tile_hm %>% save_iheatmap("/Users/ncolaian/Documents/barseq_bygroup/tiled_rare_counts_cog.html")
