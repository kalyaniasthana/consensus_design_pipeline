RUNNING test.py:
python test.py --help
usage: test.py [-h] [--fid FID] [--fidlist FIDLIST]


Give it a name YO!!

optional arguments:
  -h, --help            show this help message and exit
  --fid FID, -fid FID   Family ID. e.g.PF04398.fasta
  --fidlist FIDLIST, -fidlist FIDLIST
                        File containing list of family IDs. e.g.
                        accession_list.txt
  --strategy            select consensus design strateg 1 - with realignment, 2 - consensus with dashes 

Example:
   python test.py --fid P04532 for one family ID.
   python test.py --fidlist input_files/accession_list.txt
 
Input files needed:
   cwd/families/  needs to contain fid.fasta (P04532.fasta).

Output files:
different output files in cwd/refined_consensuses/ cwd/refined_alignments/ cwd/hmm_profiles/ cwd/hmm_emitted_sequences/ cwd/hmm_emitted_sequences_aligned/ cwd/hmm_consensuses/ cwd/dca_energy_plots/ cwd/plot_stats/ cwd/scores_refined_alignment/ cwd/scores_hmm_seqeunces/ cwd/combined_alignments/ 


For options python test.py --help
