# Code Challenge: Solve the Frequent Words Problem.

# Input: A string Text and an integer k.
# Output: All most frequent k-mers in Text.

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('dataset_2_13.txt', stringsAsFactors = F, header = F)

# split strings to characters
mystr <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
myk <- mydata[2,1] %>% pull() 

# create a place holder list
mylist <- vector(mode = "list", length = (length(mystr) - as.numeric(myk) + 1))

# put all the K-mer into mylist
for (i in 1:(length(mystr) - as.numeric(myk) + 1)) {
  mylist[[i]] <- mystr[i:(i + as.numeric(myk) - 1)]
}

# calculate freq for each k-mer and print out the most frequent k-mer
mylist %>% 
  map(~ str_c(., collapse = "")) %>% 
  unlist() %>% 
  table() %>% 
  data.frame() %>% 
  filter(Freq == max(Freq)) %>% 
  pull(1)


