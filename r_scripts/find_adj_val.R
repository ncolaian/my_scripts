# This code will find cogs that have two adjacent DA genomes. and report a matrix that contains the cog and the two associated genomes

da_tbl <- read.delim("/Users/ncolaian/Documents/v2_database/full_da_tbl_4640.txt", 
                     check.names = F, row.names = "grp_id")

#need to make 2 matrices
up_da <- da_tbl
up_da[up_da < 1] <- 0

down_da <- da_tbl
down_da[down_da > -1 | down_da < -1] <- 0

# Need to go through each cog and find adjacent values
cogs <- c() 
up_genome <- c()
down_genome <- c()

for ( rows in 1:nrow(up_da) ) {
  for ( cols in 1:( ncol(up_da)-1 ) ) {
    if ( up_da[rows,cols] == 1 & down_da[rows, (cols+1)] == -1 ) {
      cogs <- c(cogs, rownames(up_da)[rows])
      up_genome <- c(up_genome, colnames(up_da)[cols])
      down_genome <- c(down_genome, colnames(down_da)[(cols+1)])
    }
  }
}

diff_adj <- matrix(data = c(cogs, up_genome, down_genome), ncol = 3)
colnames(diff_adj)<- c("Cogs", "Up_genome", "Down_genome")
write.table(diff_adj, file="/Users/ncolaian/Documents/v2_database/opp_adj_cogs_wgenomes.txt", quote = F, sep = "\t")
