# Code Challenge: Solve the Clump Finding Problem (restated below). You will need to make sure that your algorithm is efficient enough to handle a large dataset.

# Clump Finding Problem: Find patterns forming clumps in a string.

# Input: A string Genome, and integers k, L, and t.
# Output: All distinct k-mers forming (L, t)-clumps in Genome.

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('dataset_4_5.txt', stringsAsFactors = F, header = F, fill = T)
# split strings to characters
mystr <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
k <- mydata[2,1] %>% pull() %>% str_split(., " ") %>% .[[1]] %>% .[1] %>% as.numeric()
L <- mydata[2,1] %>% pull() %>% str_split(., " ") %>% .[[1]] %>% .[2] %>% as.numeric()
t <- mydata[2,1] %>% pull() %>% str_split(., " ") %>% .[[1]] %>% .[3] %>% as.numeric()


# create a place holder list
mylist <- vector(mode = "list", length = (length(mystr) - L + 1))

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

candidates <- freq_table %>% count(patterns) %>% filter(n >= t) %>% pull(patterns)


# create a place holder list
mylist2 <- vector(mode = "list", length = length(candidates))


for (i in 1:length(candidates)) {
  # create functions to get window size for each candidate
  window_diff <- freq_table %>% filter(patterns == candidates[i]) %>% pull(rowname) %>% diff() 
  for (j in 1:length(window_diff)) {
    if (sum(window_diff[j:(j+t-2)]) <= L && !is.na(sum(window_diff[j:(j+t-2)]))) {
      mylist2[[i]] <- candidates[i]
    } else {
      next;
    }
  }
}

mylist2 # this returns nothing because there is no clump found in the dataset 

# example 2:
mystr = c("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA") %>% str_split(., "") %>% .[[1]]
k = 5 
L = 50 
t = 4

# create a place holder list
mylist <- vector(mode = "list", length = (length(mystr) - L + 1))

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

candidates <- freq_table %>% count(patterns) %>% filter(n >= t) %>% pull(patterns)


# create a place holder list
mylist2 <- vector(mode = "list", length = length(candidates))


for (i in 1:length(candidates)) {
  # create functions to get window size for each candidate
  window_diff <- freq_table %>% filter(patterns == candidates[i]) %>% pull(rowname) %>% diff() 
  for (j in 1:length(window_diff)) {
    if (sum(window_diff[j:(j+t-2)]) <= L && !is.na(sum(window_diff[j:(j+t-2)]))) {
      mylist2[[i]] <- candidates[i]
    } else {
      next;
    }
  }
}

