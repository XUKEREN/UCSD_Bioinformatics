# instructions  

# Mycobacterium tuberculosis (MTB) can persist in a latent state in humans for many years before causing disease. Latency has been found to be linked to hypoxia (lack  of oxygen) in the host. You suspect that genes that are activated in  hypoxia are regulated by a common transcription factor, so you collect  the upstream sequences for all of the MTB genes that are upregulated in  hypoxia, looking for the motif that corresponds to the binding site for  the transcription factor regulating these genes. Your biologist  colleague tells you that you should look at the 250 bp upstream region  of each gene (which have been conveniently compiled for you in a FASTA file named upstream250.txt -- right click and download this file). Your colleague also tells you that the motif is probably about 20 bp long.   


library(data.table)
library(tidyverse)
myfasta <- fread("./Finding_Hidden_Messages_in_DNA/upstream250.txt")


# We will begin by running three different software tools on the dataset provided. First, upload upstream250.txt to Consensus (Hertz and Stormo, 1999) (http://stormo.wustl.edu/consensus/html/Html/main.html).  Set the desired pattern width equal to 20 (keep all other parameters the same) and click "submit".After  the program has run, scroll to the bottom of the page and click "next".   Under "Matrix 1", you will see 19 sequences corresponding to the  substrings of the input strings having length 20 that are generated as a motif matrix.  The elements in the column to the left of these  sequences have the form XXX/YYY, where YYY represents the starting  position of each sequence in the original string of length 250.
# 
# Provide all of the starting positions of the strings of length 20.

Consensus_out <- fread("./Finding_Hidden_Messages_in_DNA/Consensus_out.txt", header = F, sep = " ")

Consensus_out  %>% separate(V4, c("junk", "start_position"), sep = "/") %>% pull(start_position) %>% noquote() 

# In order to visualize the information contained in these sequences, we will copy them into WebLogo (http://weblogo.threeplusone.com) to generate a motif logo.
# 
# Upload the image file obtained after generating this motif logo (with default parameters).  

Consensus_out %>% select(V5) %>% fwrite("Consensus_out_19_seq.txt")


# We will perform similar tasks with MEME (Bailey and Elkan, 1994) (https://meme-suite.org/meme/tools/meme).  Upload upstream250.txt, and tell MEME to find 1 motif instead of 3.  Then click on advanced options and change the minimum width to 20 and the maximum width to 20.  After submitting the process, click on "MEME html output". Notice that the motif logo has been generated under "Discovered Motifs".  Click the down arrow under "more" to see the starting positions of motifs.  To download the motif logo, click the right arrow above the logo, navigate to the "Download-logo" tab, and click "Download".
# 
# If the queue on the MEME server is too long, you can use alternate instance. 
# 
# Indicate the starting positions of the substrings of length 20 identified by MEME.

MEME_out <- fread("./Finding_Hidden_Messages_in_DNA/MEME_out.txt", header = F, sep = "\t")
MEME_out %>% pull(V4) %>% noquote()  

# Did the programs generate similar motifs?  Provide a brief (1-2 sentence) explanation.  
MEME_out %>% pull(V7) %>% noquote() -> MEME_motifs  
Consensus_out %>% pull(V5) %>% noquote() -> Consensus_motifs  
intersect(MEME_motifs, Consensus_motifs) # no overlapping  


# Although your biologist colleague told you that the motif is probably  about 20 bp long, you are skeptical, so you decide to run a motif  finding program that finds a motif over a wide range of different  lengths.
# 
# Run MEME again on upstream250.txt, but this time, use the default parameters for minimum width (6) and maximum width (50).  Note: this process may take a few minutes to run.
# 
# (a) How long is the motif produced by MEME?
#   20
# (b) Is the motif logo produced by MEME similar to the one produced before for a motif of length 20? 
#   not similar 

# (a) 2 points: The motif has length 20.
# 
# (b) 2 points: This motif is the same as the previous motif, except that it has been shifted over by a single letter.
#   
# It is better to start with shorter motifs.  As we have seen, it  takes less time to find shorter motifs, and some algorithms may even be  unable of identifying longer motifs.  Since longer motifs will contain  shorter motifs as substrings, we may be able to find longer motifs by  first finding shorter motifs and then attempting to expand them into  longer motifs.




# To evaluate the statistical significance of an identified motif, we need to ensure that a motif with the same or even larger score is unlikely to occur in a collection of "typical" DNA strings (of the same length).
# 
# How would you generate these strings?  Justify your answer.   

# There are several possible answers. Here are two:
#   
#   consider other known sequences of the same length having no motifs
# 
# randomly generate strings (ideally having the same GC-content as the species in question).
# 
# If  the motif in question has a very low probability of occurring in  randomly generated strings (or a low frequency in the known sequences),  we can conclude that it is statistically significant.


# We will now compare the different motif logos generated from varying the length of upstream regions.
# 
# Which of the motif logos that you created are similar to the motif logo generated from upstream250.txt?

# The motifs produced by upstream100, upstream500, and upstream1000 are  all similar to the motif produced by upstream250, but the motifs  produced by upstream25 does not resemble the others. (1 point for  including each of upstream100, upstream500, and upstream1000; 1 point  for not including upstream25).


