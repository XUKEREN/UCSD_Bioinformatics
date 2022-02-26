library(data.table)
library(tidyverse)

# randomly generate 500 1000-length vectors with AGCT
set.seed(1)
mylist <- list()
for (i in 1:500) {
  mylist[[i]] <- sample(c("A", "G", "C", "T"), prob = c(0.25, 0.25, 0.25, 0.25), replace = T, 1000)
}
mylist

## there are 4^9 different combinations of 9-mers

# count the number of 9-mers
mystr <- mylist[[1]]
for (i in 1:(length(mystr) - as.numeric(myk) + 1)) {
  temp[[i]] <- mystr[i:(i + as.numeric(myk) - 1)]  # store all the 9-mer from one 1000-string to a list 
}
temp %>% length() # there are 992 9-mers in a random 1000-length vector 

# the probability that the 992 is a specific 9-mer that we want to find is:
992/(4^9)

# we repeat this 500 times 
500*992/(4^9)
