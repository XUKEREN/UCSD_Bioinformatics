# We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if pi ≠ qi. For example, CGAAT and CGGAC have two mismatches. The number of mismatches between strings p and q is called the Hamming distance between these strings and is denoted HammingDistance(p, q).
# 
# Hamming Distance Problem: Compute the Hamming distance between two strings.
# 
# Input: Two strings of equal length.
# Output: The Hamming distance between these strings.

# Sample Input:
#   
# GGGCCGTTGGT
# GGACCGTTGAC
# Sample Output:
#   
# 3
library(data.table)
library(tidyverse)

HammingDistance <- function(string1, string2) {
  # create two vectors from two strings
  string1 <- string1 %>% str_split(., "") %>% .[[1]]
  string2 <- string2 %>% str_split(., "") %>% .[[1]]
  # compare two vectors 
  if (is.na((string1==string2) %>% table() %>% .["FALSE"])){
    return(0) # if the two vectors are the exact same, then return 0
  } else {
    return((string1==string2) %>% table() %>% .["FALSE"])
  }
}

# apply this function on sample input
HammingDistance("GGGCCGTTGGT" , "GGACCGTTGAC")

# apply it to the dataset 
# read dataset 
mydata = fread('dataset_9_3.txt', stringsAsFactors = F, header = F)
HammingDistance(mydata[1,1] %>% pull() , mydata[2,1] %>% pull())
