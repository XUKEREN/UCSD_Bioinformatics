# Input: An integer k, followed by a space-separated collection of strings Dna.
# Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers. (If there are multiple such strings Pattern, then you may return any one.)
# Sample Input:
#   
#   3
# AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTTCGGGACAG
# Sample Output:
#   
#   GAC

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

# get the DistanceBetweenPatternAndStrings test sample 
DistanceBetweenPatternAndStrings("AAA", c("AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"))


bases <- c("A", "C", "G", "T") 
k = 3
mykmer_list <- combn(rep(bases, k), m = k, simplify = F) 
distance_list <- mykmer_list %>% map(~ str_c(., collapse = "")) %>% map(DistanceBetweenPatternAndStrings, c("AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"))

distance_vector <- distance_list %>% unlist() 

# find the minimum value of the HammingDistance
min_value <- min(distance_vector)
# find the min index   
min_index <- which(distance_vector == min(distance_vector)) %>% min()
# find the k-mer that gave the minmum value 
min_kmer <- mykmer_list[[min_index]] 


MedianString <- function(k, Dna) {
  
  # generate all the k-mers 
  bases <- c("A", "C", "G", "T") 
  mykmer_list <- combn(rep(bases, k), m = k, simplify = F) 
  
  distance_list <- mykmer_list %>% map(~ str_c(., collapse = "")) %>% map(DistanceBetweenPatternAndStrings, Dna)
  
  distance_vector <- distance_list %>% unlist() 
  
  # find the minimum value of the HammingDistance
  min_value <- min(distance_vector)
  # find the min index   
  min_index <- which(distance_vector == min(distance_vector)) %>% min()
  # find the k-mer that gave the minmum value 
  min_kmer <- mykmer_list[[min_index]] 
  return(min_kmer)
  
}

MedianString(3, c("AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTTCGGGACAG"))


MedianString(7, c("CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC", "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"))

# modify the MedianString function to grab all kmers of all the sequences instead of listing all the potential k-mers to save time
MedianString <- function(k, Dna) {
  
  
  # a function to grab all kmers of all the sequences 
  grab_kmer <- function(x, k) {
    mystr <- x %>% str_split(., "") %>% .[[1]]
    
    # create a place holder list
    mylist <- vector(mode = "list", length = (length(mystr) - k + 1))
    
    # put all the K-mer into mylist
    for (i in 1:(length(mystr) - k + 1)) {
      mylist[[i]] <- mystr[i:(i + k - 1)]
    }
    return(mylist)
    
  }
  # apply the function to the input dna to get k-mer list 
  mykmer_list <- Dna %>% map(grab_kmer, k) %>% flatten()
  
  # calculate the distance of each k-mer with the strings 
  distance_list <- mykmer_list %>% map(~ str_c(., collapse = "")) %>% map(DistanceBetweenPatternAndStrings, Dna)
  
  distance_vector <- distance_list %>% unlist() 
  
  # find the minimum value of the HammingDistance
  min_value <- min(distance_vector)
  # find the min index   
  min_index <- which(distance_vector == min(distance_vector)) %>% min()
  # find the k-mer that gave the minmum value 
  min_kmer <- mykmer_list[[min_index]] 
  return(min_kmer)
  
}

myk <- fread("dataset_158_9.txt", fill = T) %>% .[1,1] %>% pull() %>% as.numeric()
myDna <- fread("dataset_158_9.txt", fill = T) %>% pull() %>% .[-1]
MedianString(myk, myDna)
