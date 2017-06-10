'''
Created on 2011/11/5

@author: markku
'''

from Bio import SeqIO
import random
import math

class motif_finding():
    
#initialize the motif object

    def __init__(self):
        self.seq={}
        self.seq_name={}
        self.stop_digit = {}
        self.motif_len=1 
        self.digit={}
        self.alphabet={'nt':'XACGT','aa':'ARNDCQEGHILKMFPSTWYV'}
        self.alpha_map={}
        self.type = 'nt'
        self.min_dis = 9999999 
        self.best_p = ''
        self.locat = ''
        self.local_locat_tmp = ''
        
    def load_seq(self,filename):
        print "now, I am loading your sequences ~~~ "  
        read = open(filename, 'r')
        count = 0;  
        for line in read:
            print line
            line = line.strip('\n')
            line = line.strip('')
            self.seq[count] = line
            self.seq_name[count]=''
            count = count+1
 
    def load_fasta_seq(self,filename):
        print "now, I am loading your fasta sequences ~~~ !!"
        #read seuqences from the file with 
        try:
            input_file = open(filename)
            #read sequence
            count = 0 
            test_seq = self.alphabet['aa']
            for i in range(0,4):
                test_seq.replace(self.alphabet['nt'],'')
           
                
            for seq_record in SeqIO.parse(input_file, "fasta"):
                print seq_record.id
                print repr(seq_record.seq)
                self.seq[count]=seq_record.seq
                self.seq_name[count]=seq_record.id                    
                count = count+1
            #probing sequence
            check_index = random.randint(0,count-1)
            for i in range(0,len(test_seq)):
                if (self.seq[check_index].find(test_seq[i])):
                    self.type='aa' 
                    break    
            input_file.close()
        except:
            print "reading file error"

    def set_motif_len(self,length):
        for i in range(0,length):
            #self.digit[i] = 'self.alphabet[self.type][0]'
            print i
            self.digit[i] = 1
            self.stop_digit[i] = 4
        self.motif_len = length  
        print self.motif_len 

#simulate the counter and walk through all leaves (the following functions)

    def next_leaf(self):
        #print "visit next node ...."
        for i in range(len(self.digit)-1,-1,-1):
            #print i
            if (self.digit[i] < len(self.alphabet[self.type])-1):
                #print "haha %s " % self.digit[i]
                self.digit[i] = self.digit[i]+1
                #a = a + 1, digit is actually the array
                return True
            self.digit[i] = 1
        return True    
                

    def trans_seq(self,get_digit):
        #print "transform to real AA sequence ...."
        seq = ""        
        for i in range(0,len(get_digit)):
            #print get_digit[i]
            #print self.alphabet[self.type][get_digit[i]]
            seq = seq+self.alphabet[self.type][get_digit[i]]
        #print seq
        return seq

    def cal_pair_dis(self,s,s1):
        dis = 0
        for k in range(0,self.motif_len):
                    if (s[k] != s1[k]): dis = dis + 1  
        return dis

#get total distance
    
    def total_dis(self,pat):
        total_dis = 0
        self.local_locat_tmp = ''
        #check each sequence
        self.local_locat_tmp = '['
        for i in range(0,len(self.seq)):
            #load the first motif
            mot = ''
            local_min_dis = 9999999
            cur_dis = 0
            pos = 0
            
            for j in range(0,self.motif_len):
                mot = mot + self.seq[i][j]
            local_min_dis = self.cal_pair_dis(mot,pat)
            pos = 0
             
            #check each position
            
            for j in range(self.motif_len,len(self.seq[i])):
                mot = mot +self.seq[i][j]
                mot = mot[1:]
                            
                #match the motif
               
                cur_dis = self.cal_pair_dis(mot,pat)    
                if (cur_dis < local_min_dis): 
                    local_min_dis = cur_dis
                    pos = j-self.motif_len+1
                    
            #add the minimal distance on ith sequence
            total_dis += local_min_dis
            #print "%s seq ->  %s %s" % (i, pos, local_min_dis)
            #self.local_locat_tmp += str(i)+':'+str(pos)+','
            if ( i != len(self.seq)-1):
                self.local_locat_tmp += str(pos)+','
            else:
                self.local_locat_tmp += str(pos)+''
        self.local_locat_tmp += ']'
        return total_dis
    

#Brute force median string

    def motif_finding_bfw(self):
        print "I am brute forced searching now ...."
        #    1 bestWord AAA...A 
        #    2 bestDistance 
        #   3 for each l-mer word from AAA...A to TTT...T 
        #   4     if totalDistance(word, DNA) < bestDistance
        #   5         bestDistance totalDistance(word, DNA) 
        #   6         bestWord word 
        #   7 return bestWord 
        dep = 1
        pattern = ''
        stop_pattern = self.trans_seq(self.stop_digit)

        while (dep == 1):
            
            pattern = self.trans_seq(self.digit)
            #print pattern
            temp_dis = self.total_dis(pattern)
            #print temp_dis
            if (temp_dis <= self.min_dis):
                self.min_dis = temp_dis
                self.locat=self.local_locat_tmp
                self.best_p = pattern
                #dump current result
                print "%s %s %s " % (self.best_p,self.min_dis,self.locat)
            #stopping condition
            if (pattern == stop_pattern): 
                dep = 0
            self.next_leaf()  
            #exit(); 
            #print "%s %s %s " % (self.best_p,self.min_dis,self.locat) 
        return self.best_p   
        
    def get_best_p(self):
            return self.best_p
        
    def get_min_dis(self):
            return self.min_dis
        
    def get_locat(self):
            return self.locat
        
    def get_seq(self):
            return self.seq    
    
    def reset_seq(self):
            self.seq={}
            
    def reset_all(self):
            self.stop_digit = {}
            self.motif_len=1 
            self.digit={}
            self.min_dis = 9999999 
            self.best_p = ''
            self.locat = ''
            self.local_locat_tmp = ''