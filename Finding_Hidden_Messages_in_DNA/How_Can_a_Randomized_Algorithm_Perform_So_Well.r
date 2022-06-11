# Exercise Break: Compute the probability that ten randomly selected 15-mers from the ten 600-nucleotide long strings in the Subtle Motif Problem capture at least one implanted 15-mer. (Allowable error: 0.000001)    


# 1. First you compute p1 - probability of not capturing the implanted k-mer (15-mer) in one string.

p1 = (600- 15) / (600 - 15 + 1)


# 2. Then you notice for the entire problem we have to deal with ten similar cases, i.e. you have to multiply p1 * p2... *p10, where p1 = p2 = ... = p10. So you just compute p1 to the 10th power:
pC = p1^10

# 3. Then you just compute the 'opposite' probability, i.e. the probability that from ten 600-length nucleotide string, we capture at least one implanted 15-mer! 
pAnswer = 1 - pC

# Exercise Break: Compute the probability that ten randomly selected 15-mers from ten 600-nucleotide long strings (as in the Subtle Motif Problem) capture at least two implanted 15-mers. (Allowable error: 0.000001)  

# method 1 

# 1. We have already computed p1, the probability of not capturing the 15-mer in 1st string.
p1 = (600- 15) / (600 - 15 + 1)

# 2. The opposite probability p2 is the probability of selecting correct 15-mer from 1st string:
p2 = 1 - p1

# 3. We have to consider how many there are possibilities of choosing 2 k-mers from 10 sequences - i.e. the number of combinations. For python, I did something like:

# using required libraries
library(combinat)

# generating combinations of the 
# alphabets taking 2 at a time
print ("Combination of letters two at a time")
ccounter = ombn(1:10, 2) %>% dim() %>% .[2]

# 4. To have a correct answer, multiply p2 to the power of selected k-mers, p1 to the power of 'left alone' sequences from dna, and the number of combinations:
p2^2 * p1^8 * counter


# method 2  

# In the previous exercise we calculated the probability that we will get AT LEAST one of the implanted 15-mers. So all we need to do is subtract the probability that we will get EXACTLY one 15-mer, and we will have the probability of two or more.
# 
# The probability that we will select exactly one 15-mer is the probability that we will find it on the first string and NOT on any of the others, times the probability that we will get it on the second string and NOT on the others, etc. So ten times the probability of selecting the 15-mer on the first string only.
# 
# The probability of selecting it on the first string is 1/586. The probability of NOT selecting it is thus 585/586. So for selecting it ONLY on the first string it's:
# 
# 1/586 * ((585/586) ^ 9)
# 
# Multiply by ten and subtract from the earlier result, and you're there.
pAnswer - 1/586 * ((585/586) ^ 9)*10
