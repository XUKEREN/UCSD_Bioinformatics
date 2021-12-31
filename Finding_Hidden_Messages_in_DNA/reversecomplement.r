# Reverse Complement Problem: Find the reverse complement of a DNA string.

# Input: A DNA string Pattern.
# Output: Patternrc , the reverse complement of Pattern.  

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('dataset_3_2.txt', stringsAsFactors = F, header = F)

# split strings to letters
mystr <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]

# create a place holder list
mylist <- vector(mode = "list", length = length(mystr))

# add complement 
for (i in 1:length(mystr)) {
  if (mystr[i] == "A") {
    mylist[[i]] <- "T"
  } else if (mystr[i] == "T") {
    mylist[[i]] <- "A"
  } else if (mystr[i] == "G") {
    mylist[[i]] <- "C"
  } else if (mystr[i] == "C") {
    mylist[[i]] <- "G"
  } else {
    mylist[[i]] <- "N"
  }
}

# collapse all the letters and get the reversed vector  
mylist %>% unlist() %>% rev() %>% str_c(collapse = "")

