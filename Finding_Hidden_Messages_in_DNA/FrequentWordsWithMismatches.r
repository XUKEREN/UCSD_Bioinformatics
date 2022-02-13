# One way to solve the Frequent Words with Mismatches problem is to generate all 4k k-mers Pattern, compute ApproximatePatternCount(Text, Pattern, d) for each k-mer Pattern, and then find k-mers with the maximum number of approximate occurrences. This is an inefficient approach in practice, since many of the 4k k-mers should not be considered because neither they nor their mutated versions (with up to d mismatches) appear in Text.
# 
# Instead, the following pseudocode will generalize the BetterFrequentWords function and its use of the frequency table. It uses a single map that counts the number of times a given string has an approximate match in Text. For a given k-mer substring Pattern of Text, we need to increase 1 to the count of every k-mer that has Hamming distance at most d from Pattern.  The collection of all such k-mers is called the d-neighborhood of Pattern, denoted Neighbors(Pattern, d). Check out Charging Station: Generating the Neighborhood of a String to learn how to implement Neighbors.  


# FrequentWordsWithMismatches(Text, k, d)
# Patterns ← an array of strings of length 0
# freqMap ← empty map
# n ← |Text|
#   for i ← 0 to n - k
# Pattern ← Text(i, k)
# neighborhood ← Neighbors(Pattern, d)
# for j ← 0 to |neighborhood| - 1
# neighbor ← neighborhood[j]
# if freqMap[neighbor] doesn't exist
#                 freqMap[neighbor] ← 1
#             else
#                 freqMap[neighbor] ← freqMap[neighbor] + 1
#     m ← MaxMap(freqMap)
#     for every key Pattern in freqMap
#         if freqMap[Pattern] = m
#             append Pattern to Patterns
#     return Patterns

# load the test data 
mydata <- fread("/Users/kerenxu/UCSD_Bioinformatics/Finding_Hidden_Messages_in_DNA/FrequentWordsMismatches/inputs/input_1.txt", stringsAsFactors = F, header = F)

fread("/Users/kerenxu/UCSD_Bioinformatics/Finding_Hidden_Messages_in_DNA/FrequentWordsMismatches/inputs/input_1.txt", stringsAsFactors = F)
# split strings to characters
mypattern <- mydata[1,1] %>% pull() %>% str_split(., "") %>% .[[1]]
mystring <- mydata[2,1] %>% pull() %>% str_split(., "") %>% .[[1]]
d <- mydata[3,1] %>% pull()

