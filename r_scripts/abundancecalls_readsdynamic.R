# This code will go through and answer some questions about the database mapping results
# This code creates a graph looking at the reads 

da_tbl <- read.delim("/Users/ncolaian/Documents/v2_database/full_da_tbl_4640.txt", 
                     check.names = F)

tot_reads <- read.delim("/Users/ncolaian/Documents/v2_database/reads_per_genome.txt",
                        check.names = F)

full_meta <- read.delim("/Users/ncolaian/Documents/DESeq2_test/metadata.txt")

# I want to get the top genomes 1st change all -2's and -3's to zeros
no_factors <- da_tbl[2:ncol(da_tbl)]
no_factors[no_factors <= -1] <- 0
low <- da_tbl[2:ncol(da_tbl)]
low[low != -1] <- 0

#Now need to get two vecotrs for each genome and 2 for cogs.

up_by_genome <- as.matrix(colSums(no_factors))
down_by_genome <- as.matrix(colSums(low))

#genomes
top_10_da_up <- rownames(up_by_genome)[order(up_by_genome, decreasing = T)][1:10]
bot_10_da_down <- rownames(down_by_genome)[order(down_by_genome, decreasing = F)][1:10]
#get the cogs with a DA call
agg <- rownames(down_by_genome)[down_by_genome != 0 & up_by_genome != 0]

#cogs
row.names(no_factors) <- da_tbl$grp_id
row.names(low) <- da_tbl$grp_id

up_by_cogs <- as.matrix(rowSums(no_factors))
down_by_cogs <- as.matrix(rowSums(low))

top_10_cogs <- rownames(up_by_cogs)[order(up_by_cogs, decreasing = T)][1:10]
bot_10_cogs <- rownames(down_by_cogs)[order(down_by_cogs, decreasing = F)][1:10]

agg_cogs <- rownames(down_by_cogs)[down_by_cogs != 0 & up_by_cogs != 0]

nrow(as.matrix(down_by_cogs[down_by_cogs != 0])) # This gets the number of cogs that have an abundance call

nrow(as.matrix(down_by_genome[down_by_genome != 0])) #get top genomes

# Get the total amount of reads per genome
full <- as.matrix(colSums(tot_reads[2:ncol(tot_reads)]))
full <- full+1
full <- log10(full)
  
full_reads <- merge(full_meta[,c("Id", "Fraction")], tot_reads, by.x = "Id", by.y = "SampleID")
full_reads <- aggregate(full_reads[3:ncol(full_reads)], by = list(full_reads$Fraction), FUN = sum)
full_reads <- full_reads[,2:ncol(full_reads)]+1
full_reads <- as.matrix((full_reads[1,]/full_reads[2,]))
full_reads <- t(full_reads)
full_reads[full_reads > 3] <- 3

#Try to visualize up vs down
library(ape) #Package can work with newick trees

tree <- read.tree(file = "/Users/ncolaian/Documents/v2_database/final_tree.newick")
ordered_ids <- read.delim( "/Users/ncolaian/Documents/v2_database/ordered_ids.txt", header = F)
colnames(ordered_ids) <- "ids"

# order the genome files by the tree
down_by_genome <- as.matrix(down_by_genome[match(ordered_ids$ids, rownames(down_by_genome))])
rownames(down_by_genome) <- ordered_ids$ids

up_by_genome <- as.matrix(up_by_genome[match(ordered_ids$ids, rownames(up_by_genome))])
rownames(up_by_genome) <- ordered_ids$ids

library(plotly) # Way to creat interactive heatmaps
library(iheatmapr)

# Need to make the tree into a dendogram
#tree$edge.length[tree$edge.length == 0] <- 0.00001
#tree_calib1 <- makeChronosCalib(tree, node = "root", age.min = 1, age.max = 1, interactive = F, soft.bounds = F)
#tree_try <- chronos(tree, lambda = 0, model = "relaxed", calibration = tree_calib1)
#tree_plot <- plot.phylo(tree, show.tip.label = F)

down_by_genome <- abs(down_by_genome)
down_by_genome[down_by_genome > 30 ] <- 30
up_by_genome[up_by_genome > 30 ] <- 30

da_calls_genome_heatmap <- main_heatmap(down_by_genome, name= "DA Down", colors = "Blues") %>%
  add_col_title("DA Down", side="bottom", font=list(size=12)) %>%
  add_main_heatmap(up_by_genome, name= "DA Up", colors = "Reds") %>%
  add_col_title("Differential Abundance Calls in BK vs RZ", side="top", font=list(size=18)) %>%
  add_col_title("DA Up", side="bottom", font=list(size=12)) %>%
  add_main_heatmap(full, name = "Log10(total reads)", colors = "Greens" ) %>%
  add_col_title("Log10(Total Reads)", side="bottom", font=list(size=12)) %>%
  add_main_heatmap(full_reads, name="Rz Vs Bk reads", colors = "Purples" ) %>%
  add_col_title("RZ/BK Total", side="bottom", font=list(size=12))

da_calls_genome_heatmap

da_calls_genome_heatmap %>% save_iheatmap("/Users/ncolaian/Documents/v2_database/abundance_calls_w_reads.png")
