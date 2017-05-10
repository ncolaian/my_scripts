# This script will investigate the coloniation events. Need to pass a poolcount file
# and the metadata for the samples

#Dependencies
library(ggplot2)

#Make sure to just pass in the counts to this
get_colonization_event_matrix <- function(pool_matrix, min_count) {
  pool_matrix[pool_matrix < min_count] <- 0
  pool_matrix[pool_matrix >= min_count] <- 1
  col_events <- colSums(pool_matrix)
  colon_matrix <- data.frame("Index" = names(col_events), "Col_events" = col_events)
  return(colon_matrix)
}

#The agreed upon them for figures
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
barseq_theme()
# Read in files
# Meta
meta_file <- read.csv("/Users/ncolaian/Documents/CL21_barseq_metadata.csv", header=T)
# Pool 
pool_file <- read.delim("/Users/ncolaian/Documents/Barseq_colonization/CL21_forward.poolcount",
           sep = "\t", header = T)

#This will return the colonization counts for each 
colonization_matrix <- get_colonization_event_matrix(pool_file[,6:101], 1)

#Need to merge this with the metadata ( Need Sample_type and Day Harvested )
colonization_matrix <- merge(colonization_matrix, meta_file[c("Index.1", "Day_harvest", "Sample_type")],
                             by.x = "Index", by.y = "Index.1")

# Remove the inoculums, blanks, and N/A from sample type
colonization_matrix <- colonization_matrix[!(as.character(colonization_matrix$Sample_type) == "A_LIBINOC" |
                                               as.character(colonization_matrix$Sample_type) == "B_LIBINOC" |
                                               as.character(colonization_matrix$Sample_type) == "BLANK" |
                                               as.character(colonization_matrix$Sample_type) == "N/A"
                                               ),]

#Need to order the day levels
colonization_matrix$Day_harvest <- as.character(colonization_matrix$Day_harvest)

# get the mean colonization counts for each sample type day
mean(colonization_matrix[colonization_matrix$Sample_type == "BULK_SOIL" &
                           colonization_matrix$Day_harvest == "0",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "BULK_SOIL" &
                           colonization_matrix$Day_harvest == "6",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "BULK_SOIL" &
                           colonization_matrix$Day_harvest == "12",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "BULK_SOIL" &
                           colonization_matrix$Day_harvest == "15",][,2])
#EC
mean(colonization_matrix[colonization_matrix$Sample_type == "EC" &
                           colonization_matrix$Day_harvest == "0",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "EC" &
                           colonization_matrix$Day_harvest == "6",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "EC" &
                           colonization_matrix$Day_harvest == "12",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "EC" &
                           colonization_matrix$Day_harvest == "15",][,2])
#NS
mean(colonization_matrix[colonization_matrix$Sample_type == "NEIGH_SOIL" &
                           colonization_matrix$Day_harvest == "0",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "NEIGH_SOIL" &
                           colonization_matrix$Day_harvest == "6",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "NEIGH_SOIL" &
                           colonization_matrix$Day_harvest == "12",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "NEIGH_SOIL" &
                           colonization_matrix$Day_harvest == "15",][,2])
#RHIZO
mean(colonization_matrix[colonization_matrix$Sample_type == "RHIZ" &
                      colonization_matrix$Day_harvest == "0",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "RHIZ" &
                           colonization_matrix$Day_harvest == "6",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "RHIZ" &
                           colonization_matrix$Day_harvest == "12",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "RHIZ" &
                           colonization_matrix$Day_harvest == "15",][,2])
#RR
mean(colonization_matrix[colonization_matrix$Sample_type == "RINSE_ROOT" &
                           colonization_matrix$Day_harvest == "0",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "RINSE_ROOT" &
                           colonization_matrix$Day_harvest == "6",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "RINSE_ROOT" &
                           colonization_matrix$Day_harvest == "12",][,2])
mean(colonization_matrix[colonization_matrix$Sample_type == "RINSE_ROOT" &
                           colonization_matrix$Day_harvest == "15",][,2])

#For the major colonization events
#colonization_matrix <- colonization_matrix[!(colonization_matrix$Day_harvest == "0"),]
#colonization_matrix$Day_harvest <- factor(colonization_matrix$Day_harvest, levels = c("6","12","15"))
# Part of the reordering
colonization_matrix$Day_harvest <- factor(colonization_matrix$Day_harvest, levels = c("0","6","12","15"))
#summary(colonization_matrix)

#Get the log of the colonization counts
colonization_matrix$Col_events <- colonization_matrix$Col_events +1
colonization_matrix$Col_events <- log10(colonization_matrix$Col_events)


# Create the graph
barseq_theme #load up theme

colonization_matrix$Sample_type <- ordered(colonization_matrix$Sample_type, 
                                           c("BULK_SOIL", "NEIGH_SOIL", "RHIZ", "RINSE_ROOT", "EC"))

my_plot <- ggplot(colonization_matrix, aes(Day_harvest, Col_events)) +
  geom_boxplot( outlier.colour = NA ) +
  geom_jitter()+
  facet_wrap(~Sample_type, nrow = 2)+
  ggtitle( label = "All Estimated Colonization Events in Barseq" )+
  theme( plot.title = element_text(hjust = .5) )+
  labs( x="Day Harvested", y="Log10(Colonization Events)")
my_plot
ggsave("/Users/ncolaian/Documents/Barseq_colonization/log_all_colonization.png", my_plot, width = 7)

