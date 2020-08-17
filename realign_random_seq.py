from subprocess import Popen, PIPE
from Bio import SeqIO
import random
import os
import sys

def realign(new_refined, random_seq, out_file):
    f = open(out_file, 'w+')
    p = Popen(['mafft', '--add', random_seq, '--reorder', '--keeplength', new_refined], stdout = f, stderr = PIPE)
    p.wait();stderr, stdout = p.communicate();f.flush();f.close()
    return stdout

def write_fasta_single_seq(record_seq, filename, record_id):
    with open(filename, 'w') as f:
        f.write(record_id + '\n')
        f.write(str(record_seq))

def find_random_seq(new_combined):
    for record in SeqIO.parse(new_combined, 'fasta'):
        if 'random' in record.id:
            return record

def test_func(accession):
    refined_file = os.getcwd() + '/refined_alignments/' + accession + '_refined.fasta'

    refined_sequences = list(SeqIO.parse(open(refined_file),'fasta'))
    random_seq = random.choice(refined_sequences)

    seqs = [record for record in refined_sequences if record.id != random_seq.id]
    new_refined = os.getcwd() + '/temp_files/test_refined.fasta'
    SeqIO.write(seqs, new_refined, 'fasta')

    random_seq_original = os.getcwd() + '/temp_files/random_original.fasta'
    new_combined = os.getcwd() + '/temp_files/new_combined.fasta'
    write_fasta_single_seq(str(random_seq.seq), random_seq_original, '>>random_seq_orignal')

    random_seq_no_dashes = os.getcwd() + '/temp_files/random_no_dashes.fasta'
    random_nd = list(str(random_seq.seq))
    random_nd = [i for i in random_nd if i != '-']
    random_nd = ''.join(random_nd)
    #print(random_nd)
    write_fasta_single_seq(random_nd, random_seq_no_dashes, '>>random_no_dashes')
    realign(new_refined, random_seq_no_dashes, new_combined)

    random_seq_new = find_random_seq(new_combined)
    random_seq_realigned = os.getcwd() + '/temp_files/random_realigned.fasta'
    write_fasta_single_seq(str(random_seq_new.seq), random_seq_realigned, '>>random_seq_realigned')

    print('Original: ', random_seq.seq)
    print('Realigned:', random_seq_new.seq)

    if str(random_seq.seq) == str(random_seq_new.seq):
        return True
    return False

def call_matlab(accession):
    hmm_file = os.getcwd() +  '/hmm_emitted_sequences_aligned/' + accession + '_hmmsequences_aligned.fasta'
    new_refined = os.getcwd() + '/temp_files/test_refined.fasta'
    cwd = os.getcwd()
    #os.chdir('/usr/local/MATLAB/R2016a/bin')
    cmd = './matlab -softwareopengl -nodesktop -r ' + '"cd('+"'"+cwd+'/martin_dca'+"'" + '); '+'calculate_dca_scores('+"'"+new_refined+"','"
    cmd += hmm_file+"','"+accession+"','"+cwd+"');"+'exit"'
    process = Popen(cmd,shell=True)
    process.wait()

def main():
    accession = sys.argv[1]
    flag = test_func(accession)
    print('Original == Realigned? ', flag)
    if flag is False:
        call_matlab(accession)

main()
