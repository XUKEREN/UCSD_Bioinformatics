# One way to solve the Frequent Words with Mismatches problem is to generate all 4k k-mers Pattern, compute ApproximatePatternCount(Text, Pattern, d) for each k-mer Pattern, and then find k-mers with the maximum number of approximate occurrences. This is an inefficient approach in practice, since many of the 4k k-mers should not be considered because neither they nor their mutated versions (with up to d mismatches) appear in Text.
# 
# Instead, the following pseudocode will generalize the BetterFrequentWords function and its use of the frequency table. It uses a single map that counts the number of times a given string has an approximate match in Text. For a given k-mer substring Pattern of Text, we need to increase 1 to the count of every k-mer that has Hamming distance at most d from Pattern.  The collection of all such k-mers is called the d-neighborhood of Pattern, denoted Neighbors(Pattern, d). Check out Charging Station: Generating the Neighborhood of a String to learn how to implement Neighbors.  

# FrequentWordsWithMismatches(Text, k, d)

# load the test data 
mydata <- fread("/Users/kerenxu/UCSD_Bioinformatics/Finding_Hidden_Messages_in_DNA/FrequentWordsMismatches/inputs/input_5.txt", stringsAsFactors = F)


mydata <- fread("dataset_9_9.txt", stringsAsFactors = F)

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

# create two blank lists 
# new_pattern_list1 <- list()
# new_pattern_list2 <- list()
# new_pattern_list3 <- list()

# # create the function to check the suffix pattern and append the first letter 
# AppendFirstSymbol <- function(pattern, d) {
#   
#   FirstSymbol <- pattern[1]
#   SuffixPattern <- pattern[-1]
#   
#   SuffixPattern_Neighbors <- Neighbors(SuffixPattern, d)
#   
#   for (i in 1:length(SuffixPattern_Neighbors)) {
#     
#     if (HammingDistance(SuffixPattern_Neighbors[[i]], SuffixPattern) < d) {
#       
#       nucleotide <- c("A", "C", "G", "T")
#       diff_nucleotide <- nucleotide[which(nucleotide!=FirstSymbol)] # find all the different nucleotide
#       for (j in 1:length(diff_nucleotide)) {
#         new_pattern_list1[[length(diff_nucleotide)*(i-1) + j]] <- c(diff_nucleotide[j], SuffixPattern_Neighbors[[i]])
#         
#       }
#     }
#     else if (HammingDistance(SuffixPattern_Neighbors[[i]], SuffixPattern) == d) {
#       new_pattern_list2[[i]] <- c(FirstSymbol, SuffixPattern_Neighbors[[i]])
#     } 
#     else if (HammingDistance(SuffixPattern_Neighbors[[i]], SuffixPattern) > d) {
#       new_pattern_list3[[i]] <- NULL
#     } 
#   }
#   
#   new_pattern_list12 <- c(new_pattern_list1, new_pattern_list2, new_pattern_list3, list(pattern))
#   
#   new_pattern_list12 <- new_pattern_list12[!sapply(new_pattern_list12,is.null)] %>% unique()
#   
#   return(new_pattern_list12)
#   
# }
# 

kmer_Neighborhood <- list() 

# for (i in 1:length(kmer_list)) {
#   kmer_Neighborhood[[i]] <- AppendFirstSymbol(kmer_list[[i]], d) %>% map(paste, collapse="") %>% do.call(rbind, .) %>% as.data.frame()
# }

for (i in 1:length(kmer_list)) {
  kmer_Neighborhood[[i]] <- Neighbors(kmer_list[[i]], d) %>% map(paste, collapse="") %>% do.call(rbind, .) %>% as.data.frame()
}

kmer_Neighborhood_freq_table <- kmer_Neighborhood %>% unlist() %>% data.frame() %>% table() 
which(kmer_Neighborhood_freq_table == max(kmer_Neighborhood_freq_table)) %>% names() %>% data.frame() %>% fwrite("test.results.txt")


