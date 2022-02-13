
# Our goal is to generate the d-neighborhood Neighbors(Pattern, d), the set of all k-mers whose Hamming distance from Pattern does not exceed d. As a warm up, we will first generate the 1-neigborhood of Pattern using the following pseudocode.
# 
# ImmediateNeighbors(Pattern)
# Neighborhood ← the set consisting of single string Pattern
# for i = 1 to |Pattern|
#   symbol ← i-th nucleotide of Pattern
# for each nucleotide x different from symbol
# Neighbor ← Pattern with the i-th nucleotide substituted by x
# add Neighbor to Neighborhood
# return Neighborhood

# load packages 
library(data.table)
library(tidyverse)


# initiate a Neighborhood list
Neighborhood <- list()

# create a ImmediateNeighbor function to generate 1-neigborhood of Pattern
ImmediateNeighbors <- function(pattern) {
  for (i in 1:length(pattern)) {
    symbol <- pattern[i]
    nucleotide <- c("A", "C", "G", "T")
    diff_nucleotide <- nucleotide[which(nucleotide!=symbol)] # find all the different nucleotide
    for (j in 1:length(diff_nucleotide)) {
      pattern_tmp <- pattern
      pattern[i] <- diff_nucleotide[j] # find all the different patterns than the original pattern 
      Neighborhood[[length(diff_nucleotide)*(i-1) + j ]] <- pattern # assign the pattern to the output Neighborhood list 
      pattern <- pattern_tmp
    }
  }
  Neighborhood <- c(Neighborhood, list(pattern))
  return(Neighborhood)
}

# calculate all mismatch with hamming dist 1. then compute new mismatch with ham dist 1 from prev mismatch and recurse till d.
Neighbors <- function(pattern, d) {
  Neighborhood <- ImmediateNeighbors(pattern)
  if (length(pattern) == 1) {
    
    return(print("pattern has only one letter"))
    Neighborhood <- list(NULL)
    
  } else {
    if (d == 0) {
      Neighborhood = list(pattern)
    } else if (d == 1) {
      Neighborhood
    } else if (d > 1) {
      for (i in 1:(d-1)) {
        Neighborhood <- Neighborhood %>% map(ImmediateNeighbors) %>% flatten() %>% unique()
      }
    }
    
  }
  return(Neighborhood)
}


# In the following pseudocode, we use the notation symbol • Text to denote the concatenation of a character symbol and a string Text, e.g., A•GCATG = AGCATG.
# 
# 
# Neighbors(Pattern, d)
# if d = 0
# return {Pattern}
# if |Pattern| = 1 
# return {A, C, G, T}
# Neighborhood ← an empty set
# SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
# for each string Text from SuffixNeighbors
# if HammingDistance(Suffix(Pattern), Text) < d
# for each nucleotide x
# add x • Text to Neighborhood
# else
#   add FirstSymbol(Pattern) • Text to Neighborhood
# return Neighborhood

# create a function to append the first letter later 


HammingDistance <- function(string1, string2) {
  # compare two vectors 
  if (is.na((string1==string2) %>% table() %>% .["FALSE"])){
    return(0) # if the two vectors are the exact same, then return 0
  } else {
    return((string1==string2) %>% table() %>% .["FALSE"])
  }
}

# create two blank lists 
new_pattern_list1 <- list()
new_pattern_list2 <- list()

# crreate the function to check the suffix pattern and append the first letter 
AppendFirstSymbol <- function(pattern, d) {
  
  FirstSymbol <- pattern[1]
  SuffixPattern <- pattern[-1]
  
  SuffixPattern_Neighbors <- Neighbors(SuffixPattern, d)
  
  for (i in 1:length(SuffixPattern_Neighbors)) {
    
    pattern2 <- SuffixPattern_Neighbors[[i]]
  
    if (HammingDistance(pattern2, SuffixPattern) < d) {
      
      nucleotide <- c("A", "C", "G", "T")
      diff_nucleotide <- nucleotide[which(nucleotide!=FirstSymbol)] # find all the different nucleotide
      for (j in 1:length(diff_nucleotide)) {
      new_pattern_list1[[length(diff_nucleotide)*(i-1) + j]] <- c(diff_nucleotide[j], pattern2)
        
      }
    }
    else if (HammingDistance(pattern2, SuffixPattern) == d) {
      new_pattern_list2[[i]] <- c(FirstSymbol, pattern2)
    } 
  }
  
  new_pattern_list12 <- c(new_pattern_list1, new_pattern_list2, list(pattern))
  
  new_pattern_list12 <- new_pattern_list12[!sapply(new_pattern_list12,is.null)] %>% unique()
  
  return(new_pattern_list12)
  
}


# Code Challenge: Implement Neighbors to find the d-neighborhood of a string.
# 
# Input: A string Pattern and an integer d.
# Output: The collection of strings Neighbors(Pattern, d). (You may return the strings in any order, but each line should contain only one string.)

# Sample Input:
#   
#   ACG
# 1
# Sample Output:
#   
#   CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG


# apply to the sample input 
mypattern <- "ACG" 
mypattern <- mypattern  %>% str_split(., "") %>% .[[1]]
AppendFirstSymbol(mypattern, 1)

# apply to the dataset 
mydata <- fread("dataset_3014_4.txt", stringsAsFactors = F, header = F)
# split strings to characters
mypattern <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
d <- mydata[2,1] %>% pull()
vectorList <- AppendFirstSymbol(mypattern, d) 


vectorList %>% map(paste, collapse="") %>% do.call(rbind, .) %>% as.data.frame() %>% fwrite("test.results.txt")




