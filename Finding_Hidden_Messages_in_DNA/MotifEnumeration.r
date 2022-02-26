# Implement MotifEnumeration
# 
# Input: Integers k and d, followed by a space-separated collection of strings Dna.
# Output: All (k, d)-motifs in Dna.
# MotifEnumeration(Dna, k, d)
# Patterns ← an empty set
# for each k-mer Pattern in Dna
# for each k-mer Pattern’ differing from Pattern by at most d mismatches
# if Pattern' appears in each string from Dna with at most d mismatches
#                 add Pattern' to Patterns
# remove duplicates from Patterns
# return Patterns

# Sample Input:
#   
#   3 1
# ATTTGGC TGCCTTA CGGTATC GAAAATT
# Sample Output:
#   
#   ATA ATT GTT TTT

# load libraries
library(data.table)
library(tidyverse)

# prepare the input elements 
k <- 3
d <- 1 
mystrings <- c("ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT")  

# read dataset 
k <- fread("dataset_156_8.txt", fill = T)[1,1] %>% pull() %>% as.numeric()
d <- fread("dataset_156_8.txt", fill = T)[1,2] %>% pull() %>% as.numeric()
mystrings <- fread("dataset_156_8.txt", fill = T)[2]

# initiate an output list storing all the k-mers with d mismatch for each string  
kmer_Neighborhood_output_list <- list()

# the loop function 
for (S in 1:length(mystrings)) {
  
  # split strings to characters
  mystring <- mystrings[[S]] %>% str_split(., "") %>% .[[1]]
  
  # find all the k-mers 
  
  # create a place holder list
  kmer_list <- vector(mode = "list", length = (length(mystring) - as.numeric(k) + 1))
  
  # put all the K-mer into mylist
  for (i in 1:(length(mystring) - as.numeric(k) + 1)) {
    kmer_list[[i]] <- mystring[i:(i + as.numeric(k) - 1)]
  }
  
  # find all the mismatch patterns for the k-mers 
  ###############################################
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
      Neighborhood # all the output should be in the list format 
    } else if (length(pattern) > 1) {
      if (d == 0) {
        Neighborhood = list(pattern) # all the output should be in the list format 
      } else if (d == 1) {
        Neighborhood # all the output should be in the list format 
      } else if (d > 1) {
        for (i in 1:(d-1)) {
          Neighborhood <- Neighborhood %>% map(ImmediateNeighbors) %>% flatten() %>% unique() # all the output should be in the list format 
        }
      }
      
    }
    return(Neighborhood)
  }
  
  kmer_Neighborhood <- list() 
  
  for (i in 1:length(kmer_list)) {
    kmer_Neighborhood[[i]] <- Neighbors(kmer_list[[i]], d) %>% map(paste, collapse="") %>% do.call(rbind, .) %>% as.data.frame()
  }
  
  kmer_Neighborhood_output_list[[S]] <- kmer_Neighborhood %>% unlist() %>% unique()
  
}

kmer_Neighborhood_output_list 

# intersect all the vectors in the list 
Reduce(intersect, kmer_Neighborhood_output_list) %>% data.frame() %>% fwrite("output.txt")
