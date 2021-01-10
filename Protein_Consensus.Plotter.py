        #### Protein Consensus Grapher ####

#### MODULES #####

import os
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import pandas as pd

### Input Files and Processing ###

print("PLEASE NOTE: The protein consensus grapher will work on clustalw (.clw) formatted alignment files only.") #print preliminary message
alignment_file_path=str(input("Please enter the absolute path of alignment file:")) #grab path for alignment file
name,extension=os.path.splitext(alignment_file_path) #handling file path to check for correct file extension.
if "clw" in extension:
    print("Alignment file accepted.") #progress message showing file is correct
else:
    exit("error: File extension is not recognised, please check this is a .clw file") #error message and code breaks


### Data processing from file ###

alignment = AlignIO.read(str(alignment_file_path), "clustal") #read input file using AlignIO
for protein in alignment: #iterate through every protein sequence
    sequence=protein.seq #grab sequence from protein using AligniO.seq
    for residue in sequence: #iterate over every residue within protein
        accepted_chars="ARNDBCEQZGHILKMFPSTWYV-X" #standard protein code plus - and x for alignments
        if residue in accepted_chars:  # if residue is one of the accepted characters, continue
            continue
        else:  #if residue is not recognised as accepted, exit code and instructional message to check
            exit("Error: Protein sequences have non-recognised characters. Is this a nucleotide sequence? Please check. ")  # terminate code and print error message

print("Reading file...", len(alignment), " sequences identified and accepted.") #Progress message to show proteins are read in and number entered.

### Consensus from Alignment ###

len_seqs=len(alignment[0].seq) #defining the length of our first sequence, which should be the same for all sequences.
no_seqs=len(alignment) #defining how many sequences are in the alignment
consensus_list=[] #empty list to append the most represented residues/charactersn in each column
consensus_percentages=[] #list for overall % across each column of alignment

for i in range(0, len_seqs):  # iterate through every residue from beginning to end of protein sequence.
    character_at_pos = []  # create empty list called characters at position
    for ii in range(0, no_seqs):  # iterate through each sequence in the alignment file
        char = alignment[ii].seq[i]  # defining char as each residue present (through .seq[i]), in each sequence [ii]
        character_at_pos.append(char)  # appending the char to empty list
    max_char = max(character_at_pos,key=character_at_pos.count)  # calculate max character at each position (column in alignment) stating the data as the character_at_pos list with the key being count
    if max_char == "-":  #  if the max character is equal to a gap character, this becomes 0
      max_char = 0 #alignment gaps are treated as score of 0
    consensus_list.append(max_char)  # append max char to consensus list to build a consensus sequence
    percent_consensus = (character_at_pos.count(max_char) / no_seqs) * 100  # calculate the consensus score as a %
    consensus_percentages.append(percent_consensus) #add all consensus scores to a list

#### Consensus sequence ######

consensus_seq=[] #empty list to develop the consensus sequence without gaps
for i in consensus_list: #iterate through every item in the consensus list, i, aka residue
    if i == 0: #if residue is 0 (converted gap character), skip
        continue
    else:
        consensus_seq.append(i) #appened every other residue to a new list
consensus_sequence=''.join(consensus_seq) #convert list to string
file=open("Consensus_sequence.txt","w") #open a file to write consensus string to
file.write(consensus_sequence) #write consensus sequence to file
file.close() #close file.


### Sliding Window for consensus score ###

window=int(input("enter the size of a sliding window as an odd number - 5 or 7 are recommended: ")) #user input sliding window size
if window % 2 == 0:
   exit("The window size provided is an even number \n Please provide an odd number") #error stipulating the number must be odd - centralised sliding window design
if window > len_seqs:
    exit("Window size is larger than the length of sequences in the alignment..\n Reminder: The sequence_length is " + str(len_seqs)) #error stipulating sliding window cannot exceed length of sequence

upper=(math.ceil(int(window)/2)) #upper range of sliding window, if n=7, upper section is 4
lower = (math.floor(int(window) / 2))  # lower range of sliding window, if n=7, lower section is 3
averaged=[] #empty list to compile average consensus score across sliding windows of sequence
for i in range(lower,len_seqs-lower,1): #creates a range starting from low end of window (i.e. 3) and stopping at the end-3 to prevent overlap
    chunk=consensus_percentages[i-lower:i+upper] #define windows, for example: takes 3 from one side of i and 4 from the other to create the chunk, moving by 1 each time
    averages=np.mean(chunk) #create average for 'chunk': sliding window section
    averaged.append(averages) #append to empty list which will contain average consensus score per window

### Plotting Consensus Sequence ###

fig,ax = plt.subplots(nrows=1,ncols=1) #fig,ax used to make single plot with easy customisation (incase user wants to change layout)
ax.set_xlabel("Alignment Length") #define x label
ax.set_ylabel("Consensus (%)") #define y label
ax.set_facecolor('whitesmoke') #background colour
ax.plot(averaged, '-', c='teal', linewidth=1.0) #generate line plot
plt.title("Protein Sequence Consensus Plot") #Provide title for graph
plt.savefig("Protein_Consensus_Plot.png") #save figure to user current directory

print("Protein consensus is calculated and plotted, please check for file 'protein_consensus_plot.png") #Progress message for user

### Consensus Score Tables ###

alignment_numbers=list(range(1,len_seqs+1)) #creating a range that spans the position of alignment - starting from 1 and +1 for python numbering
query = alignment[0].seq #assuming first sequence in alignment is the query sequence with a perfect match and using this for consensus scores

seqs_no_gaps=[] #creating an empty list to retain query sequence with no gaps.
for i in query: #for i, each residue, in query sequence
    if i == '-': #if residue is a gap character, continue
        continue
    else:
        seqs_no_gaps.append(i) #if residue is not a gap character, append to list to build the query sequence.


no_gaps_pos=list(range(1,len(seqs_no_gaps)+1)) #create a list of numbers that spans the length of the query sequence without gaps

temp_zip=list(zip(query,consensus_percentages)) #creating linked lists by zipping of the query sequence and consensus scores
consensus_for_zip=[] #empty list to append consensus scores to (excluding gaps)
for i in temp_zip: #iterate through every item in the temp zip
    if i[0] == '-':
        continue #ignore gap characters
    else:
        consensus_for_zip.append(i[1]) #append everything else to the list, creating consensus list with no gaps.


temp_zip2=list(zip(query,alignment_numbers)) #linking the lists of the query sequence and alignment positions

pos_in_alignment=[] #empty list to add the residue position numbers of amino acids (no gaps) from the alignment.
for i in temp_zip2: #iterate through each list couple
    if i[0] == '-':
        continue #ignore gap characters for residues
    else:
        pos_in_alignment.append(i[1]) #for non-gap characters, append their position in within the alignment to a list

complete_list=list(zip(no_gaps_pos,seqs_no_gaps,consensus_for_zip,pos_in_alignment)) #zip four lists together:
#the position of residues in the query sequence, the query sequence, the consensus scores, then the position of residues within the alignment.
df=pd.DataFrame(complete_list,columns=['Residue Number','Residue','Consensus score','Position in Alignment']) #convert zipped lists to dataframe and give column names
df.to_csv("residues_and_scores.csv") #save dataframe to csv

### Filtering of Consensus Scores ####

threshold=int(input("Enter consensus score threshold to filter conserved residues - recommended is 90%: ")) #user input threshold to filter high scoring residues out
if threshold > 100:
    exit("Selected threshold is too large, exceeds 100") #if user threshold exceeds 100, not possible and breaks code
if threshold < 0:
    exit("Selected threshold is too small, below 0.") #if user threshold below 0, not possible and breaks code

high_scoring_residues=[] #empty list to grab high scoring residues

for i in complete_list: #iterate through the 'complete list' - zipped list of all 4 lists
    if int(i[2]) >= threshold: #if consensus score for each item of complete list is equal to or greater than threshold, append to new list
        high_scoring_residues.append(i)
df1=pd.DataFrame(high_scoring_residues,columns=['Residue Number','Residue','Consensus score','Position in Alignment']) #convert the new list to df
df1.to_csv("filtered_residues.csv") #save df

print("Protein Consensus has been calculated and graphed. Check directory for following output files:"
        "\n Consensus_sequence.txt \n Protein_Consensus_Plot.png \n residues_and_scores.csv \n filtered_residues.csv")



