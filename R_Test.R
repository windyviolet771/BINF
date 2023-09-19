rm(list = ls())
library(dplyr)
exgr_test <- readRDS("D:/file/summer2/test/exgr_test.rds")

# Task 1: How many unique transcripts are there?
trans_count_ID = length(unique(exgr_test$transcript_id))

trans_count_name =length(unique(exgr_test$transcript_name))
if (trans_count_ID == trans_count_name) {
  print(paste("The number of unique transcript is:", trans_count_name))
}



# Task 2: How many unique exons are there?
exon_count = length(unique(exgr_test$exon_name))
print(paste("The number of unique exon is:", exon_count))

# find why number of exon_id is larger
new_table <- exgr_test[c("exon_id", "exon_name")]
shortened_table <- new_table %>%
  distinct(exon_id, .keep_all = TRUE)
s2 <- new_table %>%
  distinct(exon_name, .keep_all = TRUE)

# list those exons who has different ID in ChrX and ChrY
difference <- anti_join(shortened_table, s2, by = "exon_id")
head(difference)



# Task 3: what is the average length of an exon? What is the median length?
exon_mean = mean(exgr_test$width)
exon_med = median(exgr_test$width)
print(paste("The average length of an exon is:", exon_mean, "The median length is:",exon_med ))



# Task 4:Find the length of the introns between the exons.(length must be a positve number)

#Set the time count for the task
start_time <- Sys.time()
#seperate data by +/- strand, set this outside loop will reduce the running time 
data1 <- subset(exgr_test, strand == '+')
trans1 <- list() #store list

#loop over the positive strand to find the length of introns
for (i in 1:nrow(data1)) {
  #if rank is not 1, then it contains intron between it's last exons
  if (data1$rank[i] != 1) {
    #record intron's information(transcript information,rank,length)
    trans1$transcript_id[i] <- data1$transcript_id[i]
    trans1$transcript_name[i] <- data1$transcript_name[i]
    trans1$strand[i] <- '+'
    trans1$introns_rank[i] <- data1$rank[i-1]
    trans1$length[i] <- data1$start[i] -  data1$end[i-1] - 1
  } 
}

#loop over the negative strand to find the length of introns
data2 <- subset(exgr_test, strand == '-')
trans2 <- list()
for (i in 1:nrow(data2)) {
  #if rank is not 1, then it contains intron between its last exons
  if (data2$rank[i] != 1) {
    #record intron's information(transcript information,rank,length)
    trans2$transcript_id[i] <- data2$transcript_id[i]
    trans2$transcript_name[i] <- data2$transcript_name[i]
    trans2$strand[i] <- '-'
    trans2$introns_rank[i] <- data2$rank[i-1]
    trans2$length[i] <- data2$start[i-1] -  data2$end[i] - 1
  } 
}

# combine the table of negative and positive strands
table1 <- data.frame(trans1)
table2 <- data.frame(trans2)
combined_table <- rbind(table1, table2)
introns_table <- combined_table[complete.cases(combined_table), ]

# count time cost
end_time <- Sys.time()
time_cost <- end_time - start_time
print(paste("The time cost of this task is:", time_cost))

# view the head of the table to check format
head(introns_table)

# Write the result table to a tab-delimited text file
#write.table(introns_table, file = "Intron_length.txt", sep = "\t", row.names = FALSE)
