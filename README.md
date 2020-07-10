
###### I just pushed 5 new folders to the repository (PF00167, PF00902, PF04398, PF05065, PF12079). These folders have all important outputs for these 5 families. I ran the pipeline for each of these families 6 times each. 3 times using strategy 2 and 3 times using strategy 1. What these two strategies are is exlained under `Some important points`.

###### The files in these 5 folders with suffixes `_1` to `_3` are from strategy 2 executions and the ones with suffixes `_4` to `_6` are from strategy 1 executions.

### explanation for the profile matrix function
- (a profile matrix tells the probability of finding an amino acid at each position of the alignment)
- function profile_matrix(sequences_from_alignment):
	- pm = {}: keys are amino acids, including a dash. value for a key is a list. the length of this list is the length of the alignment
	(length of alignment = length of any sequence in sequencse_from_alignment)
	- intially pm[amino_acid] = [0 for 0 in range(length of alignment)]
	- for iter from 0 to length of sequences_from_alignment
		- current_sequence = sequences_from_alignment[iter]
		- for position from 0 to length of alignment
			- current_amino_acid = current_sequence[position]
			- pm[current_amino_acid][position] += 1 (increasing the frequeny of current amino acid at position
	- for amino_acid in pm:
		- values_of_amino_acid_key = pm[amino_acid]
		- for value in values_of_amino_acid_key:
			- value = value/number_of_sequences (number_of_sequences = length of sequences_from_alignment)
			
	- return pm
	
### explanation for the consensus sequence calculation function
- (finds an average sequence from all sequences in an alignment. it's an average because from each position, it takes the amino acid with the highest probability at that position)
- consensus_sequence = ''
- function consensus_sequence(pm, sequences_from_alignment):
	- for position from 0 to length of alignment:
		- current_list_of_probabilities = list()
		- for amino_acid in profile_matrix:
			- current_list_of_probabilities += [pm[amino_acid][position]]
		- max_value_for_this_position = max(current_list_of_probabilities)
		- max_indices = get a all indices in from current_list_of_probabilities which have max_value_for_this_position
		- max_index = max_indices[0]
		- if amino_acid[max_index] == '-' (if the highest probability amino acid at position is a dash '-')
			 - max_value = current_list_of_probabilities[max_index]
			 - if max_value < 0.5 (and if this highest probability is less than 0.5)
			 	- then find the second_largest_value in current_list_of_probabilities (this second largest probability should correspond to an amino acid, not a dash)
			 	- if second_largest_value == max_value:
			 		- max_index = max_indices[1]
			 	- else
			 		- max_index = index of second largest value in current_list_of_probabilities
			 - else (= if max_value > 0.5 i.e if '-' has a probability greater than 0.5)
			 	- continue
		- consensus_sequence += amino_acids[max_index] (where amino_acids = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y'])
		- Note: make sure than the keys(amino_acids) in pm are in the same order as they are in amino_acids list.

### This is what a single run through the pipeline looks like, for one family:

- remove dashes from all sequences in the alignment of the family
- check if number of sequences is less than 500, if yes then exit consensus.py
- **cd-hit clustering:** to reduce the number of sequences in the family
- check if number of sequences is less than 500, if yes then exit consnsus.py
- plot length distribution of the family
- 0th alignment step: aligning sequences generated after cd-hit step 
- Inside the while loop:
	- if any condition for breaking the loop is met
		(the breaking tags part is not necessary, can be removed or ignored)
		- save the current alignment(write_file) to **refined_alignment**
		- save the consensus sequence 
		- generate a **profile hmm** from the **refined_alignment**
		- use this **profile hmm** to **emit** N sequences (N is the number of sequences in the refined alignment)
		- align **hmm emitted sequences** with the refined alignment to generate a **combined_alignment**
		- break while loop
	- find a **profile matrix** from 'sequences' ('sequences' is the list format of sequences in write_file)
	- find a **consensus sequence** using 'sequences' and **profile matrix** (this consensus sequences does not have dashes)
	- find **bad sequences**
	- remove **bad sequences**
	- convert 'sequences' back to fasta format (temp_file)
	- remove dashes from temp_file and write output to out_file
	- align sequences in out_file and write output to write_file
	- find length of alignment  
(pydca_consensus is called from main() in consensus.py when the above while loop breaks) Inside pydca_consensus.pydca:
- split **combined_alignment** to get back the refined_alignment and an hmm emitted sequences in an alignmnet form (**hmm alignment**)
(train_test_partition functin can be ignored, we're not using it anymore. We are using the entire refined alignment to generate DCA parameters, so train_file == refined alignment and test_file == NULL. I'm going to use the term refined_alignment instead of train_file in this explanation. Looks like I should rename some variables at this point.)
- convert **hmm_sequences/hmm alignment** and **refined_alignment** to list format
- find a **profile matrix** and a consensus sequence (**hmm_consensus**) for **hmm_sequences** (this consensus sequences has dashes)
- write the above hmm_consensus to file
- realign the **refined_consensus** to the refined_alignment (refined_consensus is the consensus we got from the refined alignment and does not have dashes)
- remove the refined_consensus in it's aligned form from the refined_alignment and save it an a variable called **consensus_seq_aligned**
- write **hmm_consensus** and **consensus_seq_aligned** to **consensus_file**
- write_matlab_script writes a small matlab script
- call_matlab_script calls the above written matlab script
(Inputs to the matlab DCA program are : the refined_alignment, the hmm alignment and the consensus_file. The matlab program was not written by me, but one of shachi's collaborator's so it's more or less a black box me as I don't understand DCA that well)
- final output from the matlab program is a dca energy plot for one family, and also a plot stats file (this file has some numbers like mean, median, mode, x limits, consensus energy etc. of the dca energy plot). So the dca energy plot and the plot stats file are the two final outputs

### Some important points

- You must have noticed that the refined_consensus is generated without dashes. The refined_consensus is realigned to the refined_alignment to get back the dashes. I call this **strategy 1** of getting the refined_consensus
- There is another strategy, **strategy 2**: just get the refined_consensus with the dashes. 
- refined_consensuses from both strategies should be more or less the same right? or at least have similar dca energies? But this is not what we are observing. The refined_consensus from strategy 2 has very high dca energy for all families that I tested. Doesn't make any sense. That's why we thought that there's a big somewhere. 
- There are some useless/unused functions in consensus.py and pydca_consensus.py. I think it's time to get rid of those as well. 

###### Let me know if you want me to include explanations for each function as well.

- Are temp_files being overwritten before they are read? Multiple threads necessary? Popen helps here. 
- CD-hit input files (write.fasta) have sequences distributed across multiple lines. Does this affect SeqIO? 
- Where is consensus being calculated inside the while loop (line 394)? It's copying the mafft file (last stored from alignment(option,out_file, write_file)to refined_alignment. That is the consensus sequence? the consensus fuction is not called here.  Copy inside the if loop? Aren't we copying the same thing multiple times here. (Never mind. Didn't see the part below the break. 

