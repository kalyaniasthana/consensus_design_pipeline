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
Example:
   python test.py --fid P04532 for one family ID.
   python test.py --fidlist input_files/accession_list.txt
 
Input files needed:
   cwd/families/  needs to contain fid.fasta (P04532.fasta).

Output files:
temp_files/test_PF04398_refined.fasta
temp_files/test_PF04398_final_consensus.fasta 


Goes as far as generating a consensus sequence from THAT while loop


TODO:
1) HMM profile 2) HMM_emmit 3) DCA and everything beyond



For options python test.py --help
