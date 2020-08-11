#from subprocess import Popen,PIPE
import os
from Bio import SeqIO
import sys

def get_consensus(consensus_file):
    fasta_sequences = SeqIO.parse(open(consensus_file),'fasta')
    for record in fasta_sequences:
        cs = record.seq
    return cs

def run_test(strategy, accession):
    os.system('python test.py -fid ' + accession + ' --strategy ' + strategy)

def main(accession):
    consensus_file = os.getcwd() + '/refined_consensuses/' + accession + '_refined_consensus.fasta'
    run_test('1', accession)
    cs_1 = get_consensus(consensus_file)
    run_test('2', accession)
    cs_2 = get_consensus(consensus_file)
    return cs_1, cs_2

if __name__ == '__main__':
    result_file = os.getcwd() + '/temp_files/check_strategies.txt'
    accession = sys.argv[1]
    cs_1,cs_2 = main(accession)
    line = '{}\n{}'.format(cs_1, cs_2)
    flag = None
    if cs_1 == cs_2:
        flag = True
    else:
        flag = False
    with open(result_file, 'a') as f:
        f.write('\n' + accession + ' ' + str(flag) + '\n')
        f.write(line + '\n')
    print(line)
    print(flag)

