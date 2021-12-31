# Code Challenge: Solve the Pattern Matching Problem.

# Input: Two strings, Pattern and Genome.
# Output: A collection of space-separated integers specifying all starting positions where Pattern appears as a substring of Genome.


# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('dataset_3_5.txt', stringsAsFactors = F, header = F)

# split strings to characters
mypattern <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]

mygenome <- mydata[2,1] %>% pull() %>% str_split(., "") %>% .[[1]]

# create a place holder list
mylist <- vector(mode = "list", length = (length(mygenome) - length(mypattern) + 1))

# put all the K-mer into mylist
for (i in 1:(length(mygenome) - length(mypattern) + 1)) {
  mylist[[i]] <- mygenome[i:(i + length(mypattern) - 1)]
}

# find the matched k-mers and their index, index should minus 1
mylist %>% 
  map(~ str_c(., collapse = "")) %>% 
  unlist() %>% 
  data.frame() %>% 
  rownames_to_column() %>% 
  filter(`.` ==  str_c(mypattern, collapse = "")) %>% 
  mutate(index = as.numeric(rowname) - 1) %>% 
  pull(index) %>% str_c(., collapse = " ")

