# Code Challenge: Implement PatternCount (reproduced below).
# Input: Strings Text and Pattern.
# Output: Count(Text, Pattern).

# load packages 
library(data.table)
library(tidyverse)

# read dataset 
mydata = fread('dataset_2_6.txt', stringsAsFactors = F, header = F)

# split strings to characters
mypattern <- mydata[2,1] %>% pull() %>% str_split(., "") %>% .[[1]]

mystr = pattern <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]

# initalize the number of matched pattern = 0
npattern=0

# for loop to count the matched patterns 
for (i in 1:length(mystr)) {
        # check to see if the first element of the pattern matches any elements of the string
        if (mystr[i] == mypattern[1]) {
                # if these is match for the first element, then check the following elements
                nequal = 1
                for (j in 1:(length(mypattern) - 1)) {
                        if (mystr[i+j] == mypattern[1+j] && !is.na(mystr[i+j])) {nequal = nequal + 1}
                        else {
                                break;
                        }}
                 if (nequal == length(mypattern)) {
                         npattern = npattern + 1;   
                 } else {
                         next;     
                 }
                      
        } else {
                next;
               } 
}
print(npattern)


