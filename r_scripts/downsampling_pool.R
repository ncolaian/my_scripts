#Going to Downsample using the package caret
# To do this I need to downsample for each individual cog.I have a poolcount file
# that links the cogs to barcodes. I will then need to create a new file that has
# a downsampled poolcount file so that the barcodes with counts from each are used
library(caret)
library(tidyr)
library(mice) # waht I used to impute the data


#1) get the cogs in the pool count file
pool_file <- read.delim("/Users/ncolaian/Documents/barseq_downsample/CL21_t6_rz_bk_cogs.poolcount", header = T,
                        sep = "\t")

#2) use unique so the cogs show up once
cogs <- pool_file$geneid
cogs <- unique(cogs)

#3) get rid of barcodes that are seen 0 times in both
count_pf <- pool_file[!(rowSums(pool_file[,6:((length(colnames(pool_file)))-2)]) == 0),]

#This is a imputation check
count_pf[count_pf == 0 ] <- NA
tmp <- mice(count_pf[,6:7], m=2, maxit = 10,meth='pmm', seed=500)
densityplot(tmp)
pool_tmp <- complete(tmp, 1)
count_pf[,6:7] <- pool_tmp

#4) create a function that downsamples per cog
master_downsample <- function( pool, cog_file, sampling_times ) {
  barcodes <- c()

  # Create a matrix that will hold the counts
  downsampled_means = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(downsampled_means) <- c( "cog", colnames(pool_file)[6:((length(colnames(pool_file)))-2)] )
  
  for (scog in cog_file) {
    sub_pool <- pool[as.character(pool$geneid) == as.character(scog),]
    
    #convert to long format
    
    #need a new way to think about getting all the index columns
    sub_pool <- gather(sub_pool, index, counts, IT001:IT002)
    iden <- as.factor(sub_pool$index)
    count_pool <- sub_pool$counts
    
    #Go through and subsample
    x<-1
    #set_up the holding dataframe
    counts_from_sampling <- data.frame() #This will store the counts
    counts_from_sampling$index
    counts_from_sampling$counts

    while(x <= sampling_times) {
      x<-x+1 #make sure it keeps stepping up
      counts_from_sampling <- rbind( counts_from_sampling, get_downsample_counts( iden, count_pool ) )
    }
    if ( sum(as.numeric(counts_from_sampling$counts)) == 0 ) {
      next()
    }

    # Get the mean amount of counts
    counts_from_sampling <- aggregate(counts_from_sampling$counts, 
                                      by = list(counts_from_sampling$index), 
                                      FUN=median)
    counts_from_sampling$cog <- scog

    
    #make sure counts are as integers
    counts_from_sampling$x <- as.integer(counts_from_sampling$x)
    
    #need to put the data back into long format
    counts_from_sampling <- spread(counts_from_sampling, Group.1, x)
    downsampled_means <- rbind(downsampled_means, counts_from_sampling)
  }
  return(downsampled_means)
}
#   - Pull out the barcodes with that cog
#   - count the barcodes that have counts for each of the samples
get_downsample_counts <- function( id, cp ) {
  id <- id[!(cp == 0)]
  cp <- cp[!(cp == 0)]
  ds_pool <- downSample(cp, id)
  ds_pool <- tryCatch(aggregate(as.numeric(ds_pool$x), by=list(ds_pool$Class) , FUN = sum), error = function(cond) {
    mm <- matrix(0, ncol = 2, nrow = 1)
    mm <- as.data.frame(mm)
    return(mm)
  })
  colnames(ds_pool) <- c("index", "counts")
  return(ds_pool)
}


#   - While doing above need to create a vectors holding the barcodes that have counts

#perform the subsampling
trial <- master_downsample(count_pf, cogs, 500)

#merge with one barcode from the cog and get pool data
new_pool <- merge(trial, count_pf[,c("barcode", "rcbarcode", "scaffold", "strand", 
                                 "pos", "geneid")], by.x = "cog", by.y = "geneid")

#De-duplicate
new_pool <- new_pool[!duplicated(new_pool$cog),]

#need to re-order the columns to the original order
new_pool <- new_pool[,c( (ncol(new_pool)-4):(ncol(new_pool)), 2:(ncol(new_pool)-5))]

#print the new pool file
write.table(new_pool,sep = "\t", file = "/Users/ncolaian/Documents/barseq_downsample/CL21_indiv_rz_bk_cogs_500.poolcount",
            quote = F, row.names = F)


