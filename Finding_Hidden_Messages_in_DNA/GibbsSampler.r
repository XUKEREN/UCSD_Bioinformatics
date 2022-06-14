# load libraries  
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)

# To describe how GibbsSampler updates Motifs, we will need a slightly more advanced random number generator. Given a probability distribution (p1, …, pn), this random number generator, denoted Random(p1, …, pn), models an n-sided biased die and returns integer i with probability pi. For example, the standard six-sided fair die represents the random number generator Random(1/6, 1/6, 1/6, 1/6, 1/6, 1/6), whereas a biased die might represent the random number generator Random(0.1, 0.2, 0.3, 0.05, 0.1, 0.25). GibbsSampler further generalizes the random number generator by using the function Random(p1, …, pn) defined for any set of non-negative numbers, i.e., not necessarily satisfying the condition that the pi sum to 1. If the pi sum to some C > 0 instead, then Random(p1, …, pn) is defined as Random(p1/C, …, pn/C), where (p1/C, …, pn/C) is a probability distribution. For example, for (0.1, 0.2, 0.3) with 0.1 + 0.2 + 0.3 = 0.6,
# 
# Random(0.1, 0.2, 0.3) = Random(0.1/0.6, 0.2/0.6, 0.3/0.6) = Random(1/6, 1/3, 1/2).  

# We have previously defined the notion of a Profile-most probable k-mer in a string. We now define a Profile-randomly generated k-mer in a string Text. For each k-mer Pattern in Text, compute the probability Pr(Pattern | Profile), resulting in n = |Text| - k + 1 probabilities (p1, …, pn). These probabilities do not necessarily sum to 1, but we can still form the random number generator Random(p1, …, pn) based on them. GibbsSampler uses this random number generator to select a Profile-randomly generated k-mer at each step: if the die rolls the number i, then we define the Profile-randomly generated k-mer as the i-th k-mer in Text. While the pseudocode below repeats this procedure N times, in practice GibbsSampler depends on various stopping rules that are beyond the scope of this chapter.
# 
# GibbsSampler(Dna, kmerlength, numberofdnastrings, numberofiterations)
# motifs =  one a kmer selected uniformly at random from each string in Dna
# bestmotifs = motifs
# repeat numberofiterations times:
#   motifI = string from motifs selected uniformly at random
# profile = laplace profile of all strings in motifs except for motifI
# motifI = kmer selected using weighted random selection from ith string of dna
# if score(motifs) < score(bestmotifs)
# bestmotifs = motifs
# return bestmotifs
# STOP and Think: Note that in contrast to RandomizedMotifSearch, which always moves from higher to lower scoring Motifs, GibbsSampler may move from lower to higher scoring Motifs. Why is this reasonable?   


# Code Challenge: Implement GibbsSampler.
# 
# Input: Integers k, t, and N, followed by a space-separated collection of strings Dna.
# Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!
#   Note: The next lesson features a detailed example of GibbsSampler, so you may wish to return to this problem later.
# 
# Debug Datasets
# 
# Sample Input:
#   
#   8 5 100
# CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA
# Sample Output:
#   
#   TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG


# import the dataset - case 1
mydata <- scan("./Finding_Hidden_Messages_in_DNA/GibbsSampler/inputs/input_1.txt", character())
mystring <- mydata %>% .[-(1:3)] %>% noquote() 
k <- mydata %>% .[1] %>% as.numeric() # k-mer
t <- mydata %>% .[2] %>% as.numeric() # number of k-mer (number of the DNA string) for the profile  
N <- mydata %>% .[3] %>% as.numeric() # number of iteration 


# import the dataset - case 2
mydata <- scan("./Finding_Hidden_Messages_in_DNA/GibbsSampler/inputs/input_2.txt", character())
mystring <- mydata %>% .[-(1:3)] %>% noquote() 
k <- mydata %>% .[1] %>% as.numeric() # k-mer
t <- mydata %>% .[2] %>% as.numeric() # number of k-mer (number of the DNA string) for the profile  
N <- mydata %>% .[3] %>% as.numeric() # number of iteration 
N <- 100
Repeat <- 100


# exam dataset
mydata <- scan("./Finding_Hidden_Messages_in_DNA/dataset_163_4.txt", character())
mystring <- mydata %>% .[-(1:3)] %>% noquote() 
k <- mydata %>% .[1] %>% as.numeric() # k-mer
t <- mydata %>% .[2] %>% as.numeric() # number of k-mer (number of the DNA string) for the profile  
N <- mydata %>% .[3] %>% as.numeric() # number of iteration 
N <- 200
Repeat <- 100



# functions that will be used for GibbsSampler  


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

best_motifs_next <- profile_extractor(mystring)

# function to randomly select the k-mer for each DNA string in mystring
RandomMotifGenerator <- function(string, k) {
  
  # create a random number generator function 
  RandomNumberGenerator <- function(string, k){
    sample(1:(length(string %>% str_split(., "") %>% .[[1]]) - k + 1), 1, replace = F) # should set replace = F
  }
  
  MyRandomNumber <- RandomNumberGenerator(string, k)
  
  string %>% str_split(., "") %>% .[[1]] %>% .[MyRandomNumber:(MyRandomNumber+k-1)] %>% paste(collapse = "")
  
}



# get the probability of each kmer in string t and randomly select one kmer using weighted probability 
Profile_randomly_generated_kmer <- function(mystring) {

  # the profile from all the strings except for t
  profile <- profile_mystring_except_t %>% t()
  rownames(profile) <- c("A", "T", "C", "G")
  
  # generate all the k-mer 
  # create a place holder list
  mylist <- vector(mode = "list", length = (nchar(mystring[MyRandomStringNumber]) - k + 1))
  
  # put all the K-mer into mylist
  
  mystr <- mystring[MyRandomStringNumber] %>% str_split(., "") %>% .[[1]]
  
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
  
  # weighted random number generator  
  Profile_randomly_generated_kmer_index <- sample(1: (mylist %>% length()), 1, replace = F, prob = myprob_vector/sum(myprob_vector)) # should set replace = F
  
  # return the final selected motif 
  mylist[[Profile_randomly_generated_kmer_index]] %>% paste(collapse = "")
  
}




GibbsSampler <- function(mystring, k, t, N) {
  
  # apply the function to mystring to extract a random k-mer from each dna string.
  best_motifs <- mystring %>% map(RandomMotifGenerator, k) %>% unlist()
    
    # inner loop 
    for (j in 1:N) {
      
      MyRandomStringNumber <-  sample(1:t, 1, replace = F)
      best_motifs_except_t <- best_motifs[-MyRandomStringNumber]
      # laplace profile of all strings in motifs except for motifI
      profile_mystring_except_t <- profile_motif(best_motifs_except_t)
      
      # For each of the k-mer in the t-th sequence from the DNA, calculate Pr(k-mer | Profile).  
      # This way you will have n = |sequence|-k+1 probabilities. Use these probabilities to choose one k-mer out of n possible k-mers from the t-th sequence. 
      
      # function to weighted randomize motif  
      Profile_randomly_generated_kmer <- function(mystring) {
        
        # the profile from all the strings except for t
        profile <- profile_mystring_except_t %>% t()
        rownames(profile) <- c("A", "T", "C", "G")
        
        # generate all the k-mer 
        # create a place holder list
        mylist <- vector(mode = "list", length = (nchar(mystring[MyRandomStringNumber]) - k + 1))
        
        # put all the K-mer into mylist
        
        mystr <- mystring[MyRandomStringNumber] %>% str_split(., "") %>% .[[1]]
        
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
        
        # weighted random number generator  
        Profile_randomly_generated_kmer_index <- sample(1: (mylist %>% length()), 1, replace = F, prob = myprob_vector/sum(myprob_vector)) # should set replace = F
        
        # return the final selected motif 
        mylist[[Profile_randomly_generated_kmer_index]] %>% paste(collapse = "")
        
      }
      
      Motifs <- best_motifs
      Motifs[MyRandomStringNumber] <- Profile_randomly_generated_kmer(mystring)
      
      if (score_motif(Motifs) < score_motif(best_motifs) ) {
        best_motifs <- Motifs
      }
      
    }
  
  return(best_motifs)
  
}

# repeat_best_motifs <- list()
# repeat_best_motifs_scores <- list() 


# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)

# set up my combine function  
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

myoutput <- foreach(i=1:Repeat, .combine='comb',  .multicombine=TRUE, .init=list(list(), list()), .packages=c("data.table", "tidyverse")) %dopar% {
  repeat_best_motifs  <- GibbsSampler(mystring, k, t, N)
  repeat_best_motifs_scores <- score_motif(repeat_best_motifs)
  list(repeat_best_motifs, repeat_best_motifs_scores) # a list to hold two output
}

#stop cluster
stopCluster(cl)

repeat_best_motifs_scores <- myoutput[[2]] %>% unlist()
myoutput[[1]][which(repeat_best_motifs_scores == min(repeat_best_motifs_scores))] %>% .[[1]] %>% noquote()
