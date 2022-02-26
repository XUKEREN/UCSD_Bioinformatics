library(data.table)
library(tidyverse)

# load a motif
motif_matrix <- c(
  "TCGGGGGTTTTT",
  "CCGGTGACTTAC",
  "ACGGGGATTTTC",
  "TTGGGGACTTTT",
  "AAGGGGACTTCC",
  "TTGGGGACTTCC",
  "TCGGGGATTCAT",
  "TCGGGGATTCCT",
  "TAGGGGAACTAC",
  "TCGGGTATAACC"
)

# initiate count vector for AGCT
count_A <- rep(0, nchar(motif_matrix[1]))
count_T <- rep(0, nchar(motif_matrix[1]))
count_C <- rep(0, nchar(motif_matrix[1]))
count_G <- rep(0, nchar(motif_matrix[1]))

# j is the number of motifs in the matrix
# i is the length of each motif

for (i in 1:nchar(motif_matrix[1])) {
  for (j in 1:length(motif_matrix)) {
    
    if (motif_matrix[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "A") {
      count_A[i] = count_A[i] + 1
    } else if (motif_matrix[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "G") {
      count_G[i] = count_G[i] + 1 
    } else if (motif_matrix[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "C") {
      count_C[i] = count_C[i] + 1 
    } else if (motif_matrix[j] %>% str_split(., "") %>% .[[1]] %>% .[[i]] == "T") {
      count_T[i] = count_T[i] + 1 
    }
  }
}

motif_count <- cbind(count_A, count_T, count_C, count_G) %>% data.frame()

motif_profile <- motif_count/length(motif_matrix)

motif_score <- motif_count %>% rowSums  - motif_count %>% do.call(pmax, .)

motif_consensus <- max.col(motif_count) %>% factor(., levels = c("1", "2" , "3",   "4"), labels = c("A", "T", "C", "G"))

# the entropy of a motif matrix is defined as the sum of the entropies of its columns 
motif_profile %>% round(1) %>% mutate(
  log2count_A = ifelse(count_A == 0,0,log2(count_A)),
  log2count_T = ifelse(count_T == 0,0,log2(count_T)),  
  log2count_C = ifelse(count_C == 0,0,log2(count_C)),  
  log2count_G = ifelse(count_G == 0,0,log2(count_G)),  
  entropy = - (count_A*log2count_A + count_T*log2count_T + count_C*log2count_C + count_G*log2count_G)) %>% pull(entropy) %>% sum()
