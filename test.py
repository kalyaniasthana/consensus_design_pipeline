from subprocess import Popen, PIPE
import os
import Bio
from Bio import SeqIO,AlignIO
from ugly_strings import *
from collections import Counter, OrderedDict
import matplotlib.pyplot as plt
import numpy as np
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
        #print (x)
        return os.path.exists(x)
    #removing unwanted characters from a filename
    def refine_filename(self,ip):
        ip = str(ip, 'utf-8')
        ip = ip.strip('\n')
        ip = ip.replace('./','')
        return ip
    def list_to_file(self,list1,outfile):
        print (">>list_to_file: Write fasta without dashes.\n",outfile)
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
        print (">>Alignment:family_to_string:",fam_file,"\n")
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
        print (">>Alignment:remove_dashes_list:\n")
        #wd: list of strings with dashes
        #wod: list of strings without dashes.
        wod=[]
        seq_length=[]
        for i in wd:
            j=''.join(i).replace("-","") 
            wod+=[j];seq_length+=[len(j)]
            assert (j.isalpha()),'Error: Non-alphabet found in family'
        num_seq=len(seq_length)        
        return wod,num_seq

    def write_fasta(self,idlist,seqlist,fastafile):
        print (">>Alignment:write_fasta\n",fastafile)
        f = open(fastafile,"w+")
        #Write fasta file from ids and protein sequence list.
        assert(len(idlist)==len(seqlist)),'Idlist and sequence list do not match'
        
        for i in range(0,len(idlist)):
            id='>'+idlist[i];seq=''.join(seqlist[i])
            f.write('%s\n' % (id))
            f.write('%s\n' % (seq))
            
        f.close()
        return

    def fasta_to_mafft(self,in_file, out_file):
        print (">>Alignment:fasta_to_mafft\n",in_file,out_file)
        f=open(out_file,"w+")
        process=Popen(['mafft',in_file],stdout=f,stderr=PIPE)
        stdout, stderr = process.communicate()
        f.close()
        return stdout,stderr

    def cdhit(self,in_file, out_file):
        #cmd = 'cd-hit -i ' + in_file + ' -o ' + out_file + ' -T 1 -c 0.90'
        process=Popen(['cd-hit','-i',in_file,'-o',out_file,'-T','1','-c','0.90'],stdout=PIPE,stderr=PIPE)
        stdout, stderr = process.communicate()
        #print (stdout)
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
    
    def sequence_length_dist(self,seqlist):
        print (">>Aligment:sequence_length_dist\n")
        #Get sequence length distribution from sequencelist
        seq_length_list=[]
        for i in seqlist:
            seq_length_list+=[len(i)]
    #    print (seq_length_list)
        return seq_length_list
    def plot_length_dist(self,seq_length_list,name):
        x=np.arange(1,len(seq_length_list)+1,1)
        plt.scatter(x,seq_length_list)        
        plt.savefig(name+'.pdf')        

    def mode_of_list(self,sequence_lengths):
        print (">>Alignment:mode_of_list\n")
        n = len(sequence_lengths)
        data = Counter(sequence_lengths) 
        get_mode = dict(data) 
        mode = [k for k, v in get_mode.items() if v == max(list(data.values()))]
        if n == len(mode):
                return None
        else:
                return mode



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
    def check_break_conditions(self,num_seq,loa1,loa0,mode,count):
        print ('>>Consensus:check_break_conditions')
        #Add any number of conditions here
        #condition 1
        c1=(num_seq<500)
        #condition 2
        x = 0.1*mode[0]
        y = mode[0] - x
        c2=(loa1<y)
        #Condition3
        c3=(loa1==loa0) 
        c4=(count>10)#adding for testing
        bools=[c1,c2,c3,c4]
        print (bools)
        return bools

    def profile_matrix(self,sequences):
        print (">>Consensus:profile_matrix")
        #print (sequences)
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

    def consensus_sequence_nd(self,sequences, pm):
        consensus_seq = ''
        sequence_length = len(sequences[0])
        for i in range(sequence_length):
            l = []
            for aa in pm:
               l.append(pm[aa][i])
            max_value = max(l)
            indices = self.get_all_indices(l, max_value)
            index = indices[0]
            if self.amino_acids[index] == '-':
                if l[index] < 0.5:
                    second_largest_value = self.second_largest(l)
                    if second_largest_value == max_value:
                        index = indices[1]
                    else:
                        index = l.index(second_largest_value)
                else:
                    continue
            consensus_seq += self.amino_acids[index]

        return consensus_seq
    #finding second largest number in a list(found this on stack overflow)
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


    def consensus_sequence(self,sequences, pm):
    #find consensus sequence from sequences in list format (with dashes)
        consensus_seq = ''
        #pm = profile_matrix(sequences)
        sequence_length = len(sequences[0]) #length of any sequence
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
                        #else:
                        #        continue #if probability of occurence of dash is greater than 0.5 then skip adding an amino acid at that position

                consensus_seq += self.amino_acids[index] #append amino acid to consensus sequence

        return consensus_seq


    def iterate(self,mode,idlist,seqlist,num_seq,seq_length_list,write_file,refined_file):
        #idlist is ids of fasta sequences
        #seqlist is list of fasta sequences
        #num_seq is number of fasta sequences
        #seq_length_list is length of each fasta sequence.
        iteration =1 
        loa = 0
        count=0
        #break_tags.txt needs to be deleted manually each time?
        f_tag = open('/Users/sridharn/software/consensus_test_repo/temp_files/break_tags.txt','w+')
        print ("----------------")
        while True:
            a=Alignment()
            print("Iteration Number: " + str(iteration) + '*'*30)
            count=count+1
            print ("Count=",count)
            name_list,sequences,number_of_sequences=a.family_to_string(write_file)
            seq_length_list=a.sequence_length_dist(sequences)
            a=[]
            for i in sequences:
                a+=[str(''.join(i))]
            sequences=a
            length_of_alignment=seq_length_list[0]
            print ("Length of alignment=",length_of_alignment)
            #here loa is the length of alignment from the previous iteration
            #length_of_aligment is the length of alignment in the current iteration
            #Check for break conditions.
            breaks=self.check_break_conditions(num_seq,length_of_alignment,loa,mode,count)
            if True in breaks:
              #copy file to refined file
              p=Popen(['cp','-r',write_file,refined_file],stdin=PIPE,stdout=PIPE)
              p.communicate()
              #save consensus
              #get hmm from refined file
              #Use profile to emit N sequences.
              #align hmm sequences with refined file to generate hmm sequences.
              break
            pm = self.profile_matrix(sequences) #profile matrix
            cs = self.consensus_sequence_nd(sequences, pm) #consensus sequence at current iteration


def main():
    #All this goes into protocol.py
    ch=Check_files()
    a=Alignment()
    con=Consensus()
    #a.fasta_to_mafft('./write.fasta','./test_out')
    #Read list ofs sequences from file.
    home='/Users/sridharn/software/consensus_test_repo/'
    fam_file=home+'families/PF04398.fasta'
    idlist,seqlist,num_seq=a.family_to_string(fam_file)
    assert (len(idlist)==len(seqlist))
    assert (num_seq>=500),'Error: Less than 500 sequences in family'
    #Remove dashes
    wod,num_seq=a.remove_dashes_list(seqlist)
    #cdhit cluster
    a.write_fasta(idlist,wod,home+'temp_files/test_write.fasta')
    cdhit_out=home+'temp_files/test_cdhit.fasta'
    cdhit_in=home+'/temp_files/test_write.fasta'
    a.cdhit(cdhit_in,cdhit_out)
    cd_idlist,cd_seqlist,cd_num_seq=a.family_to_string(cdhit_out)
    print ("After clustering with CD-HIT, alignments reduced to",cd_num_seq,"from",num_seq)
    assert (cd_num_seq>=500),'Error: Less that 500 sequences after clustering with CD-HIT'
    #Plot length distribution
    seq_length_list=a.sequence_length_dist(cd_seqlist)   
    a.plot_length_dist(seq_length_list,home+'length_distributions/'+'PF04398_test')    
    #Get_mode.
    mode= a.mode_of_list(seq_length_list)
    #Zeroeth alignment:
    mafft_out=home+'temp_files/test_mafft.fasta'
    stdout,stderr=a.fasta_to_mafft(cdhit_out,mafft_out)
    mafft_idlist,mafft_seqlist,mafft_num_seq=a.family_to_string(mafft_out)
    mafft_seq_length_list=a.sequence_length_dist(mafft_out)
    #Now iterate.
    write_file=mafft_out
    refined_file= home+'temp_files/PF04398_refined.fasta'
    con.iterate(mode,mafft_idlist,mafft_seqlist,mafft_num_seq,mafft_seq_length_list,write_file,refined_file)
    


main()
