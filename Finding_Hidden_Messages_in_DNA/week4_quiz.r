# 5.
# Question 5
# Assume we are given the following strings Dna:
#   
#   ATGAGGTC
# 
# GCCCTAGA
# 
# AAATAGAT
# 
# TTGTGCTA
# 
# Then, assume that RandomizedMotifSearch begins by randomly choosing the following 3-mers Motifs of Dna:
#   
# GTC
# 
# CCC
# 
# ATA
# 
# GCT
# 
# What are the 3-mers after one iteration of RandomizedMotifSearch?  In other words, what are the 3-mers Motifs(Profile(Motifs), Dna)?  Please enter your answer as four space-separated strings.


library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)

k <- 3
t <- 4
mystring <- c("AAGCCAAA", "AATCCTGG", "GCTACTTG", "ATGTTTTG")


# for loop to iterate 1000 times 

# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)

# set up my combine function  
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

myoutput <- foreach(i=1, .combine='comb',  .multicombine=TRUE, .init=list(list(), list()), .packages=c("data.table", "tidyverse")) %dopar% {
  # for as long as the score of the constructed motifs keeps improving, which is exactly what RandomizedMotifSearch does. To implement this algorithm, you will need to randomly select the initial collection of k-mers that form the motif matrix Motifs. To do so, you will need a random number generator (denoted RandomNumber(N)) that is equally likely to return any integer from 1 to N. You might like to think about this random number generator as an unbiased N-sided die.
  

  # apply the function to mystring to extract a random k-mer from each dna string.
  best_motifs <- c("CCA", "CCT", "CTT", "TTG")
  
  # This list will be used to compare against other motif lists that we will build later in the algorithm
  
  # function to score motif 
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
  my_score_motif <- score_motif(best_motifs)
  
  # profile motif  - updated with pseudocounts.
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
    
    motif_count <- cbind(count_A + 1, count_T + 1, count_C + 1, count_G + 1)  # add pseudocounts here
    
    motif_profile <- motif_count/(length(motif) + 4) # add pseudocounts to denominator
    
    return(motif_profile)
  }
  
  # for each k-mer, extract the k-mers from all the rest DNA strings that are most_Probable   
  profile_extractor <- function(mystring) {
    motif_vector <- vector()
    
    for (a in 1:length(mystring)) {
      # create profile for the current motif_vector
      profile <- profile_motif(best_motifs) %>% t()
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
  best_motifs_next <- profile_extractor(mystring)
  my_score_motif_next <- score_motif(best_motifs_next)
  
  
  # run the while loop to check get the best motif 
  while ( my_score_motif > my_score_motif_next) {
    
    profile_extractor <- function(mystring) {
      motif_vector <- vector()
      
      for (a in 1:length(mystring)) {
        # create profile for the current motif_vector
        profile <- profile_motif(best_motifs) %>% t()
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
    
    my_score_motif <- score_motif(best_motifs)
    best_motifs_next <- profile_extractor(mystring)
    my_score_motif_next <- score_motif(best_motifs_next)
    
    best_motifs <- best_motifs_next
    best_motifs_scores <- score_motif(best_motifs)
  }
  
  list(best_motifs, best_motifs_scores) # a list to hold two output
  
}

#stop cluster
stopCluster(cl)

# find the smallest element in the list 
my_best_motifs_scores_vector <- myoutput[[2]] %>% unlist()
my_best_motifs_scores_index <- which(my_best_motifs_scores_vector == min(my_best_motifs_scores_vector))
# get the motif with the smallest index 
my_best_motifs_list <- myoutput[[1]] 
my_best_motifs_list[my_best_motifs_scores_index][[1]] %>% noquote()


# other questions  

# RandomizedMotifSearch** performs poorly when given a uniform profile matrix. → F

# it is possible for **RandomizedMotifSearch** to move from a collection of motifs with lower score to a collection of motifs with higher score. - F

# Which of the following motif-finding algorithms is guaranteed to find an optimum solution? In other words, which of the following are not heuristics? (Select all that apply.) - median string 


