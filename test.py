from subprocess import Popen, PIPE
import os
import Bio
from Bio import SeqIO,AlignIO
from ugly_strings import *
from collections import Counter, OrderedDict
class Check_files():
    '''Check all input files/dependencies/etc. are in the accessible/in the working directory.
       ugly_strings here?
    '''
    def __init__(self):
        #List input files that need to exist before running.
        #Useful to avoid hardcoded paths. 
        path=os.getcwd()+'/'
        self.cwd=path
    def common_files(self):
        cwd=self.cwd
        a,b,c,d=common_files()
        return cwd+a,cwd+b,cwd+c,cwd+d
    def specific_files(self,filename):
        cwd=self.cwd
        a,b,c,d,e,f=specific_files(filename)
        return cwd+a,cwd+b,cwd+c,cwd+d,cwd+e,cwd+f
    def pydca_strings(self,filename):
        cwd=self.cwd
        a,b,c,d,e,f,g,h,i,j=pydca_strings(filename)
        return cwd+a,cwd+b,cwd+c,cwd+d,cwd+e,cwd+f,cwd+g,cwd+h,cwd+i,cwd+j
    def fam_exist(self,accession):
        x=self.cwd+'families/'+accession+'.fasta'
        print (x)
        return os.path.exists(x)
    #removing unwanted characters from a filename
    def refine_filename(self,ip):
        ip = str(ip, 'utf-8')
        ip = ip.strip('\n')
        ip = ip.replace('./','')
        return ip
    def list_to_file(self,list1,outfile):
        print (">>list_to_file",outfile)
        with open(outfile, 'w+') as f:
            for item in list1:
                f.write("%s\n" % item)
        f.close()        
        return                

class Alignment():
    '''
    Class for performing string manipulation and alignments.
    Are all the infile and outfile variables for methods the same strings? 
    remove_dashes,cdhit,fasta_to_x,
    Alignment for one family.
    '''
    
    def __init__(self):
        #self.in_file=in_file
        #self.out_file=out_file
        self.amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y','-']
        return
    def family_to_string(self,fam_file):
        print (">>family_to_string:",fam_file,"\n")
        #Read all sequences in fasta family file and return all proteins as list of strings.
        fasta_sequences = SeqIO.parse(open(fam_file),'fasta')
        seqlist=[];idlist=[]
        for i in fasta_sequences:
            name, sequence = i.id, i.seq
            seqlist+=[list(sequence)]
            idlist+=[name]
        numseq=len(seqlist)
        print ("Read",numseq,"sequences from file",fam_file,"\n")
        return idlist,seqlist,numseq

    def remove_dashes_list(self,wd):
        print (">>remove_dashes_list:")
        #wd: list of strings with dashes
        #wod: list of strings without dashes.
        wod=[]
        seq_length=[]
        for i in wd:
            j=''.join(i).replace("-","") 
            wod+=[j];seq_length+=[len(j)]
            #assert (j.isalpha()),'Error: Non-alphabet found in family'
        print ("Remove dashes from strings")
        num_seq=len(seq_length)
        print (num_seq)        
        return wod,num_seq

    def fasta_to_clustalo(self,in_file,out_file):
        cmd = 'clustalo -i ' + in_file + ' -o ' + out_file + ' --force -v'
        process=Popen(cmd,stdout=PIPE,stderr=PIPE)
        stdout, stderr = process.communicate()
        return stdout,stderr
    def fasta_to_mafft(self,in_file, out_file):
        cmd = 'mafft ' + in_file + ' > ' + out_file
        process=Popen(['mafft',in_file],stdout=PIPE,stderr=PIPE)
        stdout, stderr = process.communicate()
        return stdout,stderr
    def fasta_to_muscle(self,in_file,out_file):
        cmd = 'muscle -in ' + in_file + ' -out ' + out_file
        process=Popen(['muscle','-in',in_file,out_file],stdout=PIPE,stderr=PIPE,shell=True)
        process.wait(2)
        stdout, stderr = process.communicate()
        print (stderr)
        return stdout,stderr
    def fasta_to_list(self,out_file):
        sequences = []
        name_list = []
        #SeqIO.parse in a function from the biopython module
        for record in Bio.SeqIO.parse(out_file, 'fasta'):
                name_list.append(record.id)
                sequences.append(str(record.seq).upper())
        return sequences, name_list

    def remove_dashes_file(self,fasta_file_from, fasta_file_to):
        with open(fasta_file_from) as fin, open(fasta_file_to, 'w') as fout:
                for line in fin:
                        if line.startswith('>'):
                                fout.write(line)
                        else:
                                fout.write(line.translate(str.maketrans('', '', '-')))
        return

    def cdhit(self,in_file, out_file):
        #cmd = 'cd-hit -i ' + in_file + ' -o ' + out_file + ' -T 1 -c 0.90'
        cmd = 'cd-hit, -i ' + in_file
        process=Popen(cmd,stdout=PIPE,stderr=PIPE)
        stdout, stderr = process.communicate()
        return stdout,stderr

    def stockholm_to_fasta(self,ifile, ofile):
        with open(ifile, 'r') as fin:
                with open(ofile, 'w') as fout:
                        sequences = Bio.SeqIO.parse(ifile, 'stockholm')
                        Bio.SeqIO.write(sequences, ofile, 'fasta')
        os.system('rm -rf ' + ifile)
        return

    def fasta_to_plain(self,accession, filename):
        alignment = AlignIO.read(open(filename, 'fasta'))
        sequences = [record.seq for record in alignment]
        plain_file = 'temp_files/' + accession + '_refined_noheader.txt'
        with open(plain_file, 'w') as f:
                for seq in sequences:
                        f.write(str(seq))
                        f.write('\n')
        return        

    def selex_to_fasta(self,in_file, out_file):
        with open(in_file) as fin, open(out_file, 'w') as fout:
                headers = []
                sequences = []
                for line in fin:
                        fout.write('>' + line[0:30].upper() + '\n')
                        fout.write(line[30: ].upper())
    
    def sequence_length_list(self,read_file):
        # returning a list of sequence lengths
        sequences, name_list = self.fasta_to_list(read_file)
        sequence_lengths = []
        #print(len(sequences))
        for seq in sequences:
                sequence_lengths.append(len(seq))
        return sequence_lengths
    #mode of a list
    def mode_of_list(self,sequence_lengths):
        n = len(sequence_lengths)
        data = Counter(sequence_lengths) 
        get_mode = dict(data) 
        mode = [k for k, v in get_mode.items() if v == max(list(data.values()))]
        if n == len(mode):
                return None
        else:
                return mode

    def alignment(self,option, in_file, out_file):
        if option == '1':
                self.fasta_to_clustalo(in_file, out_file)
        elif option == '2':
                self.fasta_to_mafft(in_file, out_file)
        elif option == '3':
                self.fasta_to_muscle(in_file, out_file)
        else:
                print('Invalid Option')
                exit()

    #realigning sequences to an existing alignment using mafft
    def realign(self,option, original_alignment, hmm_sequences, out_file):
        if option == '2':
                cwd = 'mafft --add ' + hmm_sequences + ' --reorder --keeplength ' + original_alignment + ' > ' + out_file
                #mafft --add new_sequences --reorder existing_alignment > output
                os.system(cwd)
                #op = subprocess.check_output(cwd, shell=True)
                #with open(out_file, 'w') as fin:
                #       for line in op:
                #               fin.write(line)
        else:
                print('Invalid input')
                exit()



class Consensus(object):
    '''
    Class containing various methods of finding a consensus sequence for one family.
    Should become its own module in due course.
    print (">>family_to_string:",fam_file)
    Should we pass 'sequences' as an instance for the entire Consensus class?
    '''

    def __init__(self):
        self.amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y','-']
        #self.sequences=sequences
        return
    def family_to_string(self,filename):
        print ('>>family_to_string',filename,'\n')
        #Read family from file and return all proteins as list of strings.
        #Copy to temp file and remove dashes. 
        fam=[]
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                fam+=record.id
        return fam
    def profile_matrix(self,sequences):
        '''
        Explain how you construct it. Put a reference, whatever. What is input, what is output?
        '''
        sequence_length = len(sequences[0]) #length of first sequences (length of all sequences is the same after alignment)
        profile_matrix = {} #profile matrix in dictionary format
        
        for acid in self.amino_acids:
                profile_matrix[acid] = [float(0) for i in range(sequence_length)] #initialise all entries in profile matrix to zero

        for i in range(len(sequences)):
                seq = sequences[i].upper() #convert sequence to upper case, just in case it isn't
                for j in range(len(seq)): #for each letter in the sequence
                    profile_matrix[seq[j]][j] += float(1) #increase frequency of the letter (seq[j]) at position j

        for aa in profile_matrix: #for amino acid in profile matrix
                l = profile_matrix[aa] #l i sthe list of frequencies associated with that amino acid
                for i in range(len(l)): #for position i in l
                        l[i] /= float(len(sequences)) #divide frequency at i by the length of the list l

        pm = OrderedDict([(x, profile_matrix[x]) for x in self.amino_acids])
        return pm


    def get_all_indices(self,l, value):
        return [i for i, val in enumerate(l) if val == value]
    
    def second_largest(self,numbers):
        count = 0
        m1 = m2 = float('-inf') 
        for x in numbers:
                count += 1
                if x > m2:
                        if x >= m1:
                                m1, m2 = x, m1
                        else:
                                m2 = x

        return m2 if count >= 2 else None

    def consensus1(self,sequences,pm):
        '''
        What should the variables sequences and pm contain?? 
        sequences is list read from families/pm the probability matrix
        returns one consensus sequence?
        '''
        consensus_seq=''
        sequence_length = len(sequences[0]) #length of first sequence in family.
        for i in range(sequence_length):
                l = []
                for aa in pm:
                        l.append(pm[aa][i]) #list of probabilities of amino acid 'aa' at every position
                max_value = max(l) #find maximum value in the above list
                indices = self.get_all_indices(l, max_value) #get all indices in the list which have the above maximum value
                index = indices[0] #get first index
            
                if self.amino_acids[index] == '-': #if amino acid at that index is a dash
                        if l[index] < 0.5: #if probability of occurence of dash is less than 0.5
                            second_largest_value = self.second_largest(l) #then find the second largest value
                            if second_largest_value == max_value: #if second largest and largest and largest values are equal, then get the second index from the list of max values
                                index = indices[1]
                            else:
                                index = l.index(second_largest_value) #get index of amino acid with second largest probability of occurence (after dash)
                        else:
                            continue
                consensus_seq += self.amino_acids[index]

        return consensus_seq
    def percentage_identity(self,consensus_fasta):
        A=Alignment()
        seqs, head = A.fasta_to_list(consensus_fasta)
        matches = 0
        seq_length = len(seqs[0])
        for i in range(seq_length):
                if seqs[0][i] == seqs[1][i]:
                        matches += 1
        pi = (matches*100)/seq_length
        return pi



class DCA(object):
    def __init__(self):
        return
    



def main():
    ch=Check_files()
    a=Alignment()
    con=Consensus()
    #a.fasta_to_mafft('./write.fasta','./test_out')
    #Read list of sequences from file.
    home='/Users/sridharn/software/consensus_test_repo/'
    fam_file='/Users/sridharn/software/consensus_test_repo/families/PF04398.fasta'
    idlist,seqlist,num_seq=a.family_to_string(fam_file)
    print (idlist)
    assert (len(idlist)==len(seqlist))
    assert (num_seq>=500),'Error: Less than 500 sequences in family'

    #Remove dashes
    wod,num_seq=a.remove_dashes_list(seqlist)
    
    #cdhit cluster
    ch.list_to_file(wod,home+'temp_files/test_write.fasta')            

#ß    a.cdhit(write_file,out_file)

main()
