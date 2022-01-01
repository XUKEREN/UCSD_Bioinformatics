#  How many different 9-mers form (500,3)-clumps in the E. coli genome? (In other words, do not count a 9-mer more than once.)

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread("E_coli.txt", stringsAsFactors = F, header = F, fill = T)
# split strings to characters
mystr <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
k <- 9
L <- 500
t <- 3


# create a place holder list
mylist <- vector(mode = "list", length = (length(mystr) - k + 1))

# put all the K-mer into mylist
for (i in 1:(length(mystr) - k + 1)) {
  mylist[[i]] <- mystr[i:(i + k - 1)]
}

# calculate freq for each k-mer and print out the most frequent k-mer
freq_table <- mylist %>% 
  map(~ str_c(., collapse = "")) %>% 
  unlist() %>% 
  data.frame() %>% 
  rownames_to_column() %>% 
  rename("patterns" = ".") %>% 
  mutate(rowname = as.numeric(rowname))

# obtain a list of patterns that have frequencies over t 
candidates <- freq_table %>% count(patterns) %>% filter(n >= t) %>% pull(patterns)

# optimized version of clumpfinding.r
# split to lists by unique pattern 
list_candidate_index <- freq_table %>% filter(patterns %in% candidates) %>% group_by(patterns) %>% group_split()
# get differences between consecutive indices
list_candidate_index_diff <- list_candidate_index %>% map(function(df) { df %>% pull(rowname) %>% diff()})
# create a function to find lists with window sizes smaller than L (note: the end of the k-mer should be covered by the window so should be + k - 1 )
find_small_window <- function(index) {
  df <- list_candidate_index_diff[[index]]
  for(i in 1:length(df)){
    if (df[i] + df[i+1] + k - 1 < L && !is.na(df[i+1] )) {
      return(index)
    }
  }
  
}
# apply the function find_small_window
mylist.final <- 1:length(list_candidate_index_diff) %>% map(find_small_window)
# count the number of lists that meet the criteria 
mylist.final %>% unlist() %>% data.table() %>% nrow()

