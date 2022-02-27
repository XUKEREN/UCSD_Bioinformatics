# The first potential issue with implementing MedianString is writing a function to compute the sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}. This task is achieved by the following pseudocode.

# DistanceBetweenPatternAndStrings(Pattern, Dna)
# k ← |Pattern|
#   distance ← 0
# for each string Text in Dna
# HammingDistance ← ∞
# for each k-mer Pattern’ in Text
# if HammingDistance > HammingDistance(Pattern, Pattern’)
# HammingDistance ← HammingDistance(Pattern, Pattern’)
# distance ← distance + HammingDistance
# return distance


# Sample Input:
#   
#   AAA
# TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT

# Sample Output:
#   
#   5

library(data.table)
library(tidyverse)

# HammingDistance function 
HammingDistance <- function(string1, string2) {
  # compare two vectors 
  if (is.na((string1==string2) %>% table() %>% .["FALSE"])){
    return(0) # if the two vectors are the exact same, then return 0
  } else {
    return((string1==string2) %>% table() %>% .["FALSE"])
  }
}

pattern <- "AAA"
Dna  <- c("TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT")


DistanceBetweenPatternAndStrings_summary <- function(pattern, Dna) {
  # create a list to hold the results 
  DistanceBetweenPatternAndStrings_summary_list <- list()
  
  # string split the pattern input 
  pattern <- pattern %>% str_split(., "") %>% .[[1]]
  
  
  for (S in 1:length(Dna)) {
    
    mystr <- Dna[S] %>% str_split(., "") %>% .[[1]]
    
    myk <- length(pattern)
    
    # create a place holder list
    mylist <- vector(mode = "list", length = (length(mystr) - as.numeric(myk) + 1))
    
    # put all the K-mer into mylist
    for (i in 1:(length(mystr) - as.numeric(myk) + 1)) {
      mylist[[i]] <- mystr[i:(i + as.numeric(myk) - 1)]
    }
    
    HammingDistance_summary <- vector()
    
    for (i in 1:length(mylist)) {
      HammingDistance_summary[i] <- HammingDistance(mylist[[i]], pattern)
    }
    
    HammingDistance_summary[which(HammingDistance_summary == min(HammingDistance_summary))]
    
    # find the minimum value of the HammingDistance
    min_value <- min(HammingDistance_summary)
    # find the first index that gives the minmum HammingDistance  
    min_index <- which(HammingDistance_summary == min(HammingDistance_summary)) %>% min()
    # find the k-mer that gave the minmum value 
    min_string <- mylist[[min_index]] 
    
    min_summary <- c(min_value = min_value, min_index = min_index, min_string = min_string)
    
    DistanceBetweenPatternAndStrings_summary_list[[S]] <- min_summary
    
    
    
  }
  
  return(DistanceBetweenPatternAndStrings_summary_list)
  
}


# create the DistanceBetweenPatternAndStrings function 
DistanceBetweenPatternAndStrings <- function(pattern, Dna) {
  
  DistanceBetweenPatternAndStrings_summary_list <- DistanceBetweenPatternAndStrings_summary(pattern, Dna)
  
  
  DistanceBetweenPatternAndStrings_summary_list %>%
    map(~ .[[1]]) %>% 
    unlist() %>% 
    as.numeric() %>% 
    sum()
}

# get the DistanceBetweenPatternAndStrings

DistanceBetweenPatternAndStrings("AAA", c("TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"))


# read in the input from the external dataset 
mypattern <- fread("dataset_5164_1.txt", fill = T) %>% names()
myDna <- fread("dataset_5164_1.txt", fill = T) %>% .[1,1] %>% str_split(., " ") %>% .[[1]]


DistanceBetweenPatternAndStrings(mypattern, myDna)

