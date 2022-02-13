# We now redefine the Frequent Words Problem to account for both mismatches and reverse complements. Recall that Patternrc refers to the reverse complement of Pattern.
# 
# Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a string.
# 
# Input: A DNA string Text as well as integers k and d.
# Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc) over all possible k-mers.


# Sample Input:
#   
#   ACGTTGCATGTCGCATGATGCATGAGAGCT
# 4 1
# Sample Output:
#   
#   ATGT ACAT

# load the test data 
mydata <- fread("dataset_9_10.txt", stringsAsFactors = F)

# split strings to characters
mystring <- names(mydata)[2] %>% str_split(., "") %>% .[[1]]
k <- mydata[1,1] %>% pull() 
d <- mydata[1,2] %>% pull()

# find all the k-mers 

# create a place holder list
kmer_list <- vector(mode = "list", length = (length(mystring) - as.numeric(k) + 1))

# put all the K-mer into mylist
for (i in 1:(length(mystring) - as.numeric(k) + 1)) {
  kmer_list[[i]] <- mystring[i:(i + as.numeric(k) - 1)]
}

# find all the reverse complement for the k-mers 
## create a function to find  the reverse complement for a pattern 
reverse_c <- function(mystring) {
  # create a place holder list
  reversecomplement <- vector(mode = "list", length = length(mystring))
  
  # add complement 
  for (i in 1:length(mystring)) {
    if (mystring[i] == "A") {
      reversecomplement[[i]] <- "T"
    } else if (mystring[i] == "T") {
      reversecomplement[[i]] <- "A"
    } else if (mystring[i] == "G") {
      reversecomplement[[i]] <- "C"
    } else if (mystring[i] == "C") {
      reversecomplement[[i]] <- "G"
    } else {
      reversecomplement[[i]] <- "N"
    }
  }
  # collapse all the letters and get the reversed vector  
  myreversecomplement <- reversecomplement %>% unlist() %>% rev()
  return(myreversecomplement)
  
}

# collapse all the letters and get the reversed vector  
myreversecomplement <- reversecomplement %>% unlist() %>% rev() %>% str_c(collapse = "")

# create a place holder list
kmer_list_reverse_c <- vector(mode = "list", length = (length(mystring) - as.numeric(k) + 1))

for (i in 1:(length(mystring) - as.numeric(k) + 1)) {
  kmer_list_reverse_c[[i]] <- reverse_c(kmer_list[[i]])
}

# merge all the Kmers and their reverse complement  
kmer_list <- c(kmer_list, kmer_list_reverse_c)

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


HammingDistance <- function(string1, string2) {
  # compare two vectors 
  if (is.na((string1==string2) %>% table() %>% .["FALSE"])){
    return(0) # if the two vectors are the exact same, then return 0
  } else {
    return((string1==string2) %>% table() %>% .["FALSE"])
  }
}

kmer_Neighborhood <- list() 

for (i in 1:length(kmer_list)) {
  kmer_Neighborhood[[i]] <- Neighbors(kmer_list[[i]], d) %>% map(paste, collapse="") %>% do.call(rbind, .) %>% as.data.frame()
}

kmer_Neighborhood_freq_table <- kmer_Neighborhood %>% unlist() %>% data.frame() %>% table() 
which(kmer_Neighborhood_freq_table == max(kmer_Neighborhood_freq_table)) %>% names() %>% data.frame() %>% fwrite("test.results.txt")


# try the exampple input  

#   
#   ACGTTGCATGTCGCATGATGCATGAGAGCT
# 4 1
# Sample Output:
#   
#   ATGT ACAT
# split strings to characters
mystring <- c("ACGTTGCATGTCGCATGATGCATGAGAGCT") %>% str_split(., "") %>% .[[1]]
k <- 4
d <- 1


