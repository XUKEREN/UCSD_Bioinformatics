
# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread("Salmonella_enterica.txt", stringsAsFactors = F, header = F, fill = T)

mydata = readLines("Salmonella_enterica.txt")[-1]

# Concatenate a vector of strings/character
mydata <- mydata %>% paste0(collapse = "")

# split strings to characters
mystring <- mydata %>% str_split(., "") %>% .[[1]]

# create the skew_c_g function
skew_c_g <- function(input) {
  # split the input string to a genome vector 
  genome <- input %>% str_split(., "") %>% .[[1]]
  # create an empty output vector 
  output_vector <- vector()
  # create the initial value of the output vector 
  output_vector[1] <- 0
  
  for (i in 1:length(genome)) {
    # check to see if the first element of the pattern matches any elements of the string
    if (genome[i] %in% c("A", "T")) {
      output_vector[i+1] <- output_vector[i] 
    } else if (genome[i] %in% c("C")) {
      output_vector[i+1] <- output_vector[i] - 1
    } else if (genome[i] %in% c("G")) {
      output_vector[i+1] <- output_vector[i] + 1
    }
  }
  return(output_vector)
}


# apply the function to Salmonella_enterica 
ori_location <- skew_c_g(mydata) 
which(ori_location == min(ori_location))

# ori location is 3764857  
# convert to 0-based indexing 
# create a window [A-500:A+500] 

mystring <- mystring[(3764856-500):(3764856+500)]

# find the most frequent 9-mer with 1 mismatch

k <- 9
d <- 1

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
which(kmer_Neighborhood_freq_table == max(kmer_Neighborhood_freq_table)) %>% names() %>% data.frame() %>% fwrite("DnaAboxes_search.out.txt")





