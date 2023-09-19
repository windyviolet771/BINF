rm(list = ls())
library(dplyr)
exgr_test <- readRDS("D:/file/summer2/test/exgr_test.rds")

#seperate data by the +/- strand
data1 <- subset(exgr_test, strand == '+')
data2 <- subset(exgr_test, strand == '-')

#positive strand
record1 <- list()
#loop over the positive list, leave the last case to deal individually
for (i in 1:(nrow(data1) - 1)) {
  # count the l1 region by case: at rank 1 or not 
  if (data1$rank[i] == 1){
    l1 <- 0
  } else {
    last_intro <- data1$start[i] -  data1$end[i-1] - 1
    if (last_intro < 200){ #if region is larger than 200 bp
      l1 <- data1$start[i] - ceiling(last_intro/2) 
    } else {
      l1 <- data1$start[i] - 100
    }
  }
  
  #count the l2 and u1 region, also see if region is larger than 200 bp
  if (data1$width[i] < 200){
    count <- ceiling(data1$width[i]/2)
    l2 <- data1$start[i] + count
    u1 <- l2 + 1
  } else {
    l2 <- data1$start[i] + 100
    u1 <- data1$end[i] - 100
  }
  
  #count the u2 regions by case: if it's the last tank insides its transcript
  if (data1$rank[i+1] == 1) {
    u2 <- 0
  } else {
    next_intro <- data1$start[i+1] -  data1$end[i] - 1
    if (next_intro < 200){
      u2 <- data1$start[i+1] - ceiling(next_intro/2) - 1
    } else {
      u2 <- data1$end[i] + 100
    }
  }
  #record l1,l2,u1,u2 regions and it's position
  record1$transcript_id[i] <- data1$transcript_id[i]
  record1$transcript_name[i] <- data1$transcript_name[i]
  record1$exon[i] <- data1$exon_name[i]
  record1$start[i] <- data1$start[i]
  record1$end[i] <- data1$end[i]
  record1$width[i] <- data1$width[i]
  record1$strand[i] <- "+"
  record1$l1[i] <- l1
  record1$l2[i] <- l2
  record1$u1[i] <- u1
  record1$u2[i] <- u2
}

# edge case(last row of the data)
# count the l2 and u2 regions at edge case
if (data1$width[nrow(data1)] < 200){
  count <- ceiling(data1$width[nrow(data1)]/2)
  l2 <- data1$start[nrow(data1)] + count
  u1 <- l2 + 1
} else {
  l2 <- data1$start[nrow(data1)] + 100
  u1 <- data1$end[nrow(data1)] - 100
}

# count the l1 and u2 regions at edge case
if (data1$rank[nrow(data1)] == 1){# the edge case is a individual transcript
  l1 <- 0
  u2 <- 0
} else {# the edge case isn't a single transcript
  
  last_intro <- data1$start[nrow(data1)] -  data1$end[nrow(data1)-1] - 1
  if (last_intro < 200){
    l1 <- data1$start[nrow(data1)] - ceiling(last_intro/2) 
  } else {
    l1 <- data1$start[nrow(data1)] - 100
  }
  u2 <- 0
}
#record l1,l2,u1,u2 regions at edge case
record1$transcript_id[nrow(data1)] <- data1$transcript_id[nrow(data1)]
record1$transcript_name[nrow(data1)] <- data1$transcript_name[nrow(data1)]
record1$exon[nrow(data1)] <- data1$exon_name[nrow(data1)]
record1$start[nrow(data1)] <- data1$start[nrow(data1)]
record1$end[nrow(data1)] <- data1$end[nrow(data1)]
record1$width[nrow(data1)] <- data1$width[nrow(data1)]
record1$strand[nrow(data1)] <- "+"
record1$l1[nrow(data1)] <- l1
record1$l2[nrow(data1)] <- l2
record1$u1[nrow(data1)] <- u1
record1$u2[nrow(data1)] <- u2



#negative strand
record2 <- list()
#loop over the positive list, leave the last case to deal individually
for (i in 1:(nrow(data2) - 1)) {
  # count the l1 region by case: at rank 1 or not 
  if (data2$rank[i] == 1){
    u2 <- 0
  } else {
    last_intro <- data2$start[i-1] -  data2$end[i] - 1
    if (last_intro < 200){
      u2 <- data2$end[i] + ceiling(last_intro/2) 
    } else {
      u2 <- data2$end[i] + 100
    }
  }
  
  #count the l2 and u1 region, also see if region is larger than 200 bp
  if (data2$width[i] < 200){
    count <- ceiling(data2$width[i]/2)
    l2 <- data2$start[i] + count 
    u1 <- l2 + 1
  } else {
    u1 <- data2$start[i] + 100
    l2 <- data2$end[i] - 100
  }
  
  #count the u2 regions by case: if it's the last tank insides its transcript
  if (data2$rank[i+1] == 1) {
    l1 <- 0
  } else {
    next_intro <- data2$start[i] - data2$end[i+1] - 1
    if (next_intro < 200){
      l1 <- data2$end[i+1] + ceiling(next_intro/2) + 1
    } else {
      l1 <- data2$start[i] - 100
    }
  }
  #record l1,l2,u1,u2 regions and it's position
  record2$transcript_id[i] <- data2$transcript_id[i]
  record2$transcript_name[i] <- data2$transcript_name[i]
  record2$exon[i] <- data2$exon_name[i]
  record2$start[i] <- data2$start[i]
  record2$end[i] <- data2$end[i]
  record2$width[i] <- data2$width[i]
  record2$strand[i] <- "-"
  record2$l1[i] <- l1
  record2$l2[i] <- l2
  record2$u1[i] <- u1
  record2$u2[i] <- u2
}

# edge case(last row of the data)
# count the l2 and u2 regions at edge case
if (data2$width[nrow(data2)] < 200){
  count <- ceiling(data2$width[nrow(data2)]/2)
  l2 <- data2$start[nrow(data2)] + count 
  u1 <- l2 + 1
} else {
  u1 <- data2$start[nrow(data2)] + 100
  l2 <- data2$end[nrow(data2)] - 100
}
# count the l1 and u2 regions at edge case
if (data2$rank[nrow(data2)] == 1){# the edge case is a individual transcript
  l1 <- 0
  u2 <- 0
} else {# the edge case isn't a single transcript
  last_intro <- data2$start[nrow(data2)-1] -  data2$end[nrow(data2)] - 1
  if (last_intro < 200){
    u2 <- data2$end[i] + ceiling(last_intro/2) 
  } else {
    u2 <- data2$end[i] + 100
  }
  l1 <- 0
}
#record l1,l2,u1,u2 regions at edge case
record2$transcript_id[nrow(data2)] <- data2$transcript_id[nrow(data2)]
record2$transcript_name[nrow(data2)] <- data2$transcript_name[nrow(data2)]
record2$exon[nrow(data2)] <- data2$exon_name[nrow(data2)]
record2$start[nrow(data2)] <- data2$start[nrow(data2)]
record2$end[nrow(data2)] <- data2$end[nrow(data2)]
record2$width[nrow(data2)] <- data2$width[nrow(data2)]
record2$strand[nrow(data2)] <- "-"
record2$l1[nrow(data2)] <- l1
record2$l2[nrow(data2)] <- l2
record2$u1[nrow(data2)] <- u1
record2$u2[nrow(data2)] <- u2

#combine files separated by +/- strand
region <- data.frame(record1)
region2 <- data.frame(record2)
combined_table <- rbind(region, region2)

#write the output file
Bouns_table <- combined_table[complete.cases(combined_table), ]
#write.table(Bouns_table, file = "Bouns_table.txt", sep = "\t", row.names = FALSE)
