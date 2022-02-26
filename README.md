# UCSD_Bioinformatics
this repo holds my scripts for the UCSD Bioinformatics Specialization on Coursera

Solving the Frequent Words Problem to find k-mers that appear several times within a window  
- `FrequentWords.r`  
- `PatternCount.r`  
- `reversecomplement.r`  
- `patternmatch.r`  
- `Apply.patternmatch.r`  
- `ClumpFinding.r`  
- `ApplyClumpFinding.r` 

Solving the Minimum Skew Problem now provides us with an approximate location of ori   
- `Skew_C_G.r`  
- `locating_ori.r`

Hamming Distance Problem  
- `HammingDistance.r`  

Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string  
- `approx_pattern_matching.r`  
- `count_approx_pattern_matching.r`  

Implement Neighbors to find the d-neighborhood of a string
- `Neighborhood.r`  

Solve the Frequent Words with Mismatches Problem  
- `FrequentWordsWithMismatches.r`   

Frequent Words with Mismatches and Reverse Complements Problem  
- `FrequentWordsWithMismatchesComplements.r`

Thermotoga petrophila DnaAboxes finder 
- `DnaAboxes_search.r`  

Bibliography Notes  
- [Sedgewick and Flajolet, 2013 gave an overview of various approaches for computing the probabilities of patterns in a string.](https://aofa.cs.princeton.edu/home/)


Which DNA patterns play the role of molecular clocks  
- What is the expected number of occurrences of a 9-mer in 500 random DNA strings, each of length 1000? Assume that the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25). `9-mer.in500randomDNAstrings.1000length.r`   

Motif finding  
- `MotifEnumeration.r`  

Subtle Motif Problem  
- The Subtle Motif Problem refers to implanting a 15-mer with four random mutations in ten randomly generated 600 nucleotide-long strings (the typical length of many upstream regulatory regions). The instance of the Subtle Motif Problem that we will use has the implanted 15-mer AAAAAAAAGGGGGGG and is given below.  [dataset](/Finding_Hidden_Messages_in_DNA/subtle_motif_dataset.txt)


Scoring Motifs  



