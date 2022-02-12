# Let's follow the 5' → 3' direction of DNA and walk along the chromosome from ter to ori (along a reverse half-strand), then continue on from ori to ter (along a forward half-strand). In our previous discussion, we saw that the skew is decreasing along the reverse half-strand and increasing along the forward half-strand. Thus, the skew should achieve a minimum at the position where the reverse half-strand ends and the forward half-strand begins, which is exactly the location of ori!
# 
# We have just developed an insight for a new algorithm for locating ori: it should be found where the skew attains a minimum.
# 
# Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
# 
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('dataset_7_10.txt', stringsAsFactors = F, header = F)

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

# apply it to mydata 
skew_c_g(mydata)

# return index of the smallest value in the vector 
a <- skew_c_g(mydata)
which.min(a)
# results needs 1-based indexing 
which(a == min(a))-1

