# Since we don't know the location of ori in a circular genome, let's linearize it (i.e., select an arbitrary position and pretend that the genome begins here), resulting in a linear string Genome. We define Skewi(Genome) as the difference between the total number of occurrences of G and the total number of occurrences of C in the first i nucleotides of Genome. The skew diagram is defined by plotting Skewi (Genome) (as i ranges from 0 to |Genome|), where Skew0 (Genome) is set equal to zero. The figure below shows a skew diagram for the DNA string CATGGGCATCGGCCATACGCC.
# 
# Note that we can compute Skewi+1(Genome) from Skewi(Genome) according to the nucleotide in position i of Genome. If this nucleotide is G, then Skewi+1(Genome) = Skewi(Genome) + 1; if this nucleotide is C, then Skewi+1(Genome)= Skewi(Genome) – 1; otherwise, Skewi+1(Genome) = Skewi(Genome).
# Sample Input:
#   CATGGGCATCGGCCATACGCC
# 
# Sample Output:
#   0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

# load libraries
library(data.table)
library(tidyverse)

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

# test the function using the sample input
skew_c_g("CATGGGCATCGGCCATACGCC")

# apply it to "GAGCCACCGCGATA"
skew_c_g("GAGCCACCGCGATA")



skew_c_g("GCATACACTTCCCAGTAGGTACTG")
