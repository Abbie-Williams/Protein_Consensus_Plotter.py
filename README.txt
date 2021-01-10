-------------------------------------------------------------------------------------------------
|				PROTEIN CONSENSUS PLOTTER			        	|
|											        |
|	AUTHOR: ABBIE WILLIAMS 								        |
|											        |
|	NOTE: THIS SCRIPT WAS DEVELOPED FOR ASSESSMENT WITHIN THE DATA SCIENCE AND MACHINE      |	
|	      LEARNING FOR THE BIOSCIENCES 2020 MODULE, AS PART OF THE SWBIO DTP PROGRAM. 	|
|										 	        |	
-------------------------------------------------------------------------------------------------


			I. PROTEIN CONSENSUS PLOTTER 
			II. KEY CONCEPTS AND ASSUMPTIONS
			III. EXPECTED WORKFLOW
			IV. DEPENDENCIES AND RUNNING THE SCRIPT  
			V. FILE INPUT 
			VI. FILE OUTPUTS
			VII. EXAMPLE FILES 
			VIII. ACKNOWLEDGEMENTS / REFERENCES 
			IX. CONTACT INFORMATION 


I. PROTEIN CONSENSUS PLOTTER 
-------------------------------------------------------------------------------------------------
The protein consensus plotter is a basic tool developed to calculate and visualise the 
conservation of amino acid residues across a multiple protein sequence alignment. The theory 
behind this is that as protein function is conserved amongst a group of related proteins, 
specific amino acid(s) determining the protein structure will also be conserved, such as 
domains, motifs, or catalytic sites. The protein consensus plotter can be used to identify well 
conserved residue(s) for several purposes: preliminary analysis of proteins of an unknown
function, analysing the conservation of residues of interest i.e. frequently mutated residues, 
or comprehensive analysis of the conservation of entire protein sequence i.e. amongst protein 
families. 

II. KEY CONCEPTS AND ASSUMPTIONS
------------------------------------------------------------------------------------------------

Prior to running, a basic understanding of the protein consensus plotter is required. The 
protein consensus plotter will analyse a protein alignment file to calculate consensus. 
Consensus is defined as 'agreement per residue': the level of conservation of the most 
represented residue for each position across the alignment sequence. Hence, consensus is 
effectively a measure of amino acid conservation and is converted to a percentage by:

	(Frequency of most represented character* / Number of sequences in alignment ) x 100

*calculated for each residue of alignment, across all sequences. 

In the instance of gap characters, "-', these are treated as a value of 0 and are plotted. Graph
plotting uses sliding windows (n=user defined) to create rolling averages which smooth the plot 
line. When considering particular residues and their consensus scores, their position and score 
can be inferred from the graph and cross referenced with the output file 'consensus_scores.csv' 
for more accurate analysis. The consensus plot should be used as a visual aid and the 
consensus_score.csv file should always be checked for more information. Importantly, a key 
assumption of the protein consensus plotter is that the first sequence of the alignment file 
is the protein query sequence. This assumption is based on the intended workflow (see III.
EXPECTED WORKFLOW for explanation). If this is not the case, the alignment file must be 
reformatted and/or realigned to meet this assumption.   

III. EXPECTED WORKFLOW
-------------------------------------------------------------------------------------------------

I. Retrieve the protein sequence of interest (FASTA)

II. Run BLASTP for protein query*. If filtering of hits are required, this should be done at 
this stage. A manual check of whether the first hit matches the protein sequence should also be
carried out. If this condition is not met as the query protein is not present in the selected
database, the protein sequence (FASTA) can be added as the first entry of the multifasta file
downloaded from BLAST. All desired hits to be used in the consensus analysis should be 
downloaded as the complete sequences in a multifasta file. 

III. Align the protein sequences downloaded from BLAST using a multiple sequence aligner (MSA).
The MSA can be selected based on user preference, and ran locally or online**. Importantly, 
the alignment file must be CLUSTALW format (.clw) for the script to run.

IV. Run the Protein Consensus Plotter script, using the alignment file as an input when prompted. 

V. Check consensus plot for runs of conserved residues (with a high consensus score). Refer back
to output files listing consensus scores to identify the exact consensus score, residue type, 
and residue position in the sequence and in the alignment. 

*available here: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins
**available here: https://www.ebi.ac.uk/Tools/msa/

IV. DEPENDENCIES AND RUNNING THE SCRIPT
-------------------------------------------------------------------------------------------------
The Protein Consensus Plotter was developed and tested with Python ver.3.7 using Pycharm 
ver.2020.2.3 and requires the following modules to run:

os
math
numpy
pandas
Biopython
Matlotib

The Protein Consensus Plotter is largely guided by user input. Simply, the script can be run by 
entering the following in a terminal:

	"python path_to_script_location/Protein_Consensus_Plotter.py"

Input files are not specified in the command, but are then prompted by the user. There are three
user input prompts: the absolute file path of the input file, the size of the sliding window to
be used, and the threshold value for filtering residues based on consensus score. All user
prompts are signaled with obvious, clear messages. 

V. FILE INPUT
-------------------------------------------------------------------------------------------------
A single  MSA alignment file in CLUSTALW format (.clw) only. See folder 'example files'.  

VI. FILE OUTPUTS
-------------------------------------------------------------------------------------------------
'Consensus_sequence.txt': The overall consensus sequence generated from the alignment (contains
the most common residues across entire alignment).
'Protein_consensus_plot.png': Graph of sequence consensus (%) across the alignment length.
'Residues_and_scores.csv': Dataframe containing the residue number (within query sequence),
Residue (amino acid), Consensus score (%) and Position in alignment (residue number in alignment
sequence).
'Filtered_residues.csv': Dataframe containing the residue number (within query sequence),
Residue (amino acid), Consensus score (%) and Position in alignment (residue number in alignment
sequence) for all residues exceeding the specified threshold.  

VII. EXAMPLE FILES
-------------------------------------------------------------------------------------------------
"alignment_input_file.clw": MSA file used for input. This is based on a cathepsin (protease) 
protein in Fasciola hepatica* which was processed in the expected workflow (see III. EXPECTED 
WORKFLOW).Cathespin protein sequence was used in a BLASTP search with default settings. 100 
sequences were then aligned with using MUSCLE with default parameters, generating the example 
input files. 

The following are all examples of the output file generated from using the input file:
"Protein_Consensus_Plot.png","Consensus_sequence.txt","filtered_residues.csv", and
"residues_and_scores.csv".


*available here: https://www.ncbi.nlm.nih.gov/protein/AAA29137.1?from=1&to=326&report=fasta


VIII. ACKNOWLEDGEMENTS / REFERENCES 
-------------------------------------------------------------------------------------------------
Pycharm: JetBrains, 2017. Pycharm. JetBrains. Available at: https://www.jetbrains.com/pycharm/
Python:Python Software Foundation. Python Language Reference, version 3.7. Available at 
http://www.python.org

os: os - miscellaneous operating system interfaces. Manual available here: 
https://docs.python.org/3/library/os.html

Math: math - Mathematical functions. Manual available here: 
https://docs.python.org/3/library/math.html

numpy: Charles R. Harris, K. Jarrod Millman, Stéfan J. van der Walt, Ralf Gommers, 
Pauli Virtanen, David Cournapeau, Eric Wieser, Julian Taylor, Sebastian Berg, Nathaniel J. Smith,
 Robert Kern, Matti Picus, Stephan Hoyer, Marten H. van Kerkwijk, Matthew Brett, Allan Haldane, 
Jaime Fernández del Río, Mark Wiebe, Pearu Peterson, Pierre Gérard-Marchant, Kevin Sheppard, 
Tyler Reddy, Warren Weckesser, Hameer Abbasi, Christoph Gohlke & Travis E. Oliphant. Array
 programming with NumPy, Nature, 585, 357–362 (2020)

pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th
Python in Science Conference, 51-56 (2010)

Biopython: Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck 
T, Kauff F, Wilczynski B and de Hoon MJL (2009) Biopython: freely available Python tools for 
computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423

Matplotlib: J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & 
Engineering, vol. 9, no. 3, pp. 90-95, 2007.


IX. CONTACT INFORMATION 
-------------------------------------------------------------------------------------------------
For queries or problems, please contact: hj20164@bristol.ac.uk
