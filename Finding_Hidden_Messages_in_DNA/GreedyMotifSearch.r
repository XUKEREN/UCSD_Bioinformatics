library(data.table)
library(tidyverse)

# Code Challenge: Implement GreedyMotifSearch.
# 
# Input: Integers k and t, followed by a space-separated collection of strings Dna.
# Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t). If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
# Note: If you are not satisfied with the performance of GreedyMotifSearch — even if you implemented it correctly — then wait until we discuss this algorithm in the next lesson.
# 
# Debug Datasets
# 
# Sample Input:
#   
#   3 5
# GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG
# Sample Output:
#   
#   CAG CAG CAA CAA CAA

# GreedyMotifSearch(Dna, k, t)

# great post explaining the Pseudo Code
library(data.table)
library(tidyverse)

# import the dataset 
mystring <- fread("dataset_159_5.txt", fill = T) %>% .[2] %>% as.vector() %>% as.character()
k <- fread("dataset_159_5.txt", fill = T) %>% .[1,1] %>% as.numeric() # k-mer
t <- fread("dataset_159_5.txt", fill = T) %>% .[1,2] %>% as.numeric() # number of k-mer (number of the DNA string) for the profile 

# first find out the first k-mer for each DNA string in mystring
BestMotifs_function <- function(string, k) {
  string %>% str_split(., "") %>% .[[1]] %>% .[1:k] %>% paste(collapse = "")
}
# apply the function to mystring to extract the first k-mer from each dna string.
best_motifs_int <- mystring %>% map(BestMotifs_function, k) %>% unlist()

# This list will be used to compare against other motif lists that we will build later in the algorithm.  

# score motif 
score_motif <- function(motif) {
  # initiate count vector for AGCT
  count_A <- rep(0, nchar(motif[1]))
  count_T <- rep(0, nchar(motif[1]))
  count_C <- rep(0, nchar(motif[1]))
  count_G <- rep(0, nchar(motif[1]))
  
  # j is the number of motifs in the matrix
  # i is the length of each motif
  
  for (i in 1:nchar(motif[1])) {
    for (j in 1:length(motif)) {
      
      if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "A") {
        count_A[i] = count_A[i] + 1
      } else if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "G") {
        count_G[i] = count_G[i] + 1 
      } else if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "C") {
        count_C[i] = count_C[i] + 1 
      } else if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "T") {
        count_T[i] = count_T[i] + 1 
      }
    }
  }
  
  motif_count <- cbind(count_A, count_T, count_C, count_G) %>% data.frame()
  
  motif_profile <- motif_count/length(motif)
  
  motif_score <- motif_count %>% rowSums  - motif_count %>% do.call(pmax, .) 
  
  motif_score_sum <- sum(motif_score)
  
  return(motif_score_sum)
}

# get the score from the profile comprising the first k-mer from each dna string 
score_motif_int <- score_motif(best_motifs_int)

# profile motif 
profile_motif <- function(motif) {
  # initiate count vector for AGCT
  count_A <- rep(0, nchar(motif[1]))
  count_T <- rep(0, nchar(motif[1]))
  count_C <- rep(0, nchar(motif[1]))
  count_G <- rep(0, nchar(motif[1]))
  
  # j is the number of motifs in the matrix
  # i is the length of each motif
  
  for (i in 1:nchar(motif[1])) {
    for (j in 1:length(motif)) {
      
      if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "A") {
        count_A[i] = count_A[i] + 1
      } else if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "G") {
        count_G[i] = count_G[i] + 1 
      } else if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "C") {
        count_C[i] = count_C[i] + 1 
      } else if (motif[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "T") {
        count_T[i] = count_T[i] + 1 
      }
    }
  }
  
  motif_count <- cbind(count_A, count_T, count_C, count_G) 
  
  motif_profile <- motif_count/length(motif)  
  
  return(motif_profile)
}

# create a function to find the profile for the first k-mer in the first string 
# init a vector to hold all the profile_most_probable motif 

# find all the k-mers for the first Dna string
first_DNA_string <- mystring[1] %>% str_split(., "") %>% .[[1]] 

# create a place holder list
first_DNA_string_kmer <- vector(mode = "list", length = (length(first_DNA_string) - k + 1))

# put all the K-mer into mylist
for (i in 1:(length(first_DNA_string) - k + 1)) {
  first_DNA_string_kmer[[i]] <- first_DNA_string[i:(i + k - 1)]
}

# get all the kmer from the first DNA string
first_DNA_string_kmer

# for each k-mer, extract the k-mers from all the rest DNA strings that are most_Probable 
profile_extractor <- function(kmer) {
  motif_vector <- vector()
  motif_vector[1] <- kmer %>% paste(collapse = "")
  
  for (a in 2:length(mystring)) {
    # create profile for the current motif_vector
    profile <- profile_motif(motif_vector) %>% t()
    rownames(profile) <- c("A", "T", "C", "G")
    
    # find out the profile_most_probable motif 
    
    # find out the profile_most_probable motif from the next Dna string
    mystr <- mystring[a] %>% str_split(., "") %>% .[[1]]
    
    # generate all the k-mer 
    
    # create a place holder list
    mylist <- vector(mode = "list", length = (length(mystr) - k + 1))
    
    # put all the K-mer into mylist
    for (i in 1:(length(mystr) - k + 1)) {
      mylist[[i]] <- mystr[i:(i + k - 1)]
    }
    
    
    Profile_most_Probable <- function(x) {
      myprob <- prod(
        profile["T",][which(x == "T")] %>% prod(),
        profile["A",][which(x == "A")] %>% prod(),
        profile["G",][which(x == "G")] %>% prod(),
        profile["C",][which(x == "C")] %>% prod()
      )
      
      return(myprob)
    }
    
    
    myprob_vector <- mylist %>% map(Profile_most_Probable) %>% unlist() 
    
    motif_vector[a] <- mylist[which(myprob_vector == max(myprob_vector))] %>% .[[1]] %>% paste(collapse = "")
    
    
  }
  
  return(motif_vector)
}

# get all the motif profiles and calculate the score 
final_motif_vector <- first_DNA_string_kmer %>% map(profile_extractor) 
final_score_vector <- final_motif_vector %>% map(score_motif) %>% unlist()

# check the minimum score to see if its smaller than score_motif_int
score_motif_int
final_score_vector[which(final_score_vector == min(final_score_vector)) %>% min() ] 

if (score_motif_int <= final_score_vector[which(final_score_vector == min(final_score_vector)) %>% min() ] ) {
  print(best_motifs_int %>% noquote()) 
} else {
  print(final_motif_vector[[which(final_score_vector == min(final_score_vector)) %>% min() ]] %>% noquote())
}
