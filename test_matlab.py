from subprocess import Popen, PIPE
import os
import sys

home = os.getcwd()
#print(home)
#sys.exit()
process = Popen(['whereis', 'MATLAB'], stdout = PIPE, stderr = PIPE)
stdout, stderr = process.communicate();process.wait()
path = str(stdout, 'utf-8').split(': ')[1].strip('\n')
os.chdir(path + '/R2016a/bin')
# process = Popen(['./matlab', '-softwareopengl', '-nodesktop', '-r', '"cd(])
fasta_train = home+'/refined_alignments/PF04398_refined.fasta'
fasta_test = home+'/hmm_emitted_sequences_aligned/PF04398_hmmsequences_aligned.fasta'
# ./matlab -softwareopengl -nodesktop -r "cd('/media/Data/consensus_project/martin_dca'); calculate_dca_scores('/media/Data/consensus_project/refined_alignments/PF04398_refined.fasta','/media/Data/consensus_project/hmm_emitted_sequences_aligned/PF04398_hmmsequences_aligned.fasta', 'PF04398', '/media/Data/consensus_project');exit"
accession = 'PF04398'
cmd = './matlab -softwareopengl -nodesktop -r ' + '"cd('+"'"+ home+'/martin_dca'+"'" + '); '+'calculate_dca_scores('+"'"+fasta_train+"','"
cmd += fasta_test+"','"+accession+"','"+home+"');"+'exit"'

#process = Popen('./matlab -softwareopengl -nodesktop -r "cd('+home+'martin_dca/); calculate_dca_scores('+fasta_train+','+fasta_test+','+accession+','+home');exit"', shell=True) 
process = Popen(cmd,shell=True)
process.wait()



