# Our goal now is to modify our previous algorithm for the Frequent Words Problem in order to find DnaA boxes by identifying frequent k-mers, possibly with mismatches. Given strings Text and Pattern as well as an integer d, we define Countd(Text, Pattern) as the total number of occurrences of Pattern in Text with at most d mismatches. For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. Note that two of these occurrences overlap.
# 
# Exercise Break: Compute Count2(AACAAGCTGATAAACATTTAAAGAG, AAAAA).

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
mypattern <- "AAAAA"
mystring <- "AACAAGCTGATAAACATTTAAAGAG"
d <- 1
# split strings to vectors 
mypattern <- mypattern %>% str_split(., "") %>% .[[1]]
mystring <- mystring %>% str_split(., "") %>% .[[1]]
approx_pattern_matching(mypattern, mystring, d) %>% length()

# apply the function to the new strings 
# read dataset 
mypattern <- "AAAAA"
mystring <- "AACAAGCTGATAAACATTTAAAGAG"
d <- 2
# split strings to vectors 
mypattern <- mypattern %>% str_split(., "") %>% .[[1]]
mystring <- mystring %>% str_split(., "") %>% .[[1]]
approx_pattern_matching(mypattern, mystring, d) %>% length()

# create a new count function
count <- function(pattern, string, d){
  approx_pattern_matching(pattern, string, d) %>% length()
}

count(mypattern, mystring, d)

# try the sample input 
mypattern <- "GAGG"
mystring <- "TTTAGAGCCTTCAGAGG"
d <- 2
mypattern <- mypattern %>% str_split(., "") %>% .[[1]]
mystring <- mystring %>% str_split(., "") %>% .[[1]]
count(mypattern, mystring, d)

# apply it to the dataset 
# read dataset 
mydata = fread('dataset_9_6.txt', stringsAsFactors = F, header = F)
# split strings to characters
mypattern <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
mystring <- mydata[2,1] %>% pull() %>% str_split(., "") %>% .[[1]]
d <- mydata[3,1] %>% pull()
count(mypattern, mystring, d)

