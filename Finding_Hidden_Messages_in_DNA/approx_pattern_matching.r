# We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear with slight variations leads to the following generalization of the Pattern Matching Problem.
# 
# Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
# 
# Input: Strings Pattern and Text along with an integer d.
# Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
# Code Challenge: Solve the Approximate Pattern Matching Problem.


# Sample Input:
#   
# ATTCTGGA
# CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
# 3
# Sample Output:
#   
# 6 7 26 27

# load packages 
library(data.table)
library(tidyverse)

approx_pattern_matching <- function(pattern, string, d) {
  # create an empty vector for the output 
  output_vector <- vector()  
  # HammingDistance function 
  HammingDistance <- function(string1, string2) {
    # compare two vectors 
    if (is.na((string1==string2) %>% table() %>% .["FALSE"])){
      return(0)
    } else {
      return((string1==string2) %>% table() %>% .["FALSE"])
    }
    
  }
  
  # loop compare the pattern and a subset of the string 
  for (i in 1:(length(string) - length(pattern) + 1) ) {
    if (HammingDistance(pattern, string[i: (i + length(pattern) - 1)]) <= d) {
      output_vector[i] <- i
    } else {
      output_vector[i] <- NA
    }
  }
  
  # return 0-based indexing 
  return(output_vector[which(!is.na(output_vector))] - 1)
}


# test the function on sample input  
mypattern <- "ATTCTGGA"
mystring <- "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
d <- 3
# split strings to vectors 
mypattern <- mypattern %>% str_split(., "") %>% .[[1]]
mystring <- mystring %>% str_split(., "") %>% .[[1]]
approx_pattern_matching(mypattern, mystring, d)

# apply the function to the dataset 
# read dataset 
mydata = fread('dataset_9_4.txt', stringsAsFactors = F, header = F)
# split strings to characters
mypattern <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
mystring <- mydata[2,1] %>% pull() %>% str_split(., "") %>% .[[1]]
d <- mydata[3,1] %>% pull()
approx_pattern_matching(mypattern, mystring, d) %>% as.list() %>% fwrite("results.txt", sep = " ")



