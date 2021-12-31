# Exercise Break: Return a space-separated list of starting positions (in increasing order) where CTTGATCAT appears as a substring in the Vibrio cholerae genome.

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('Vibrio_cholerae.txt', stringsAsFactors = F, header = F)

# split strings to characters
mypattern <- c("CTTGATCAT") %>% str_split(., "") %>% .[[1]]

mygenome <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]

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

