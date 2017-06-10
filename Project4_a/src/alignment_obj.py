'''
Created on 2011/12/1

author: markku

ID: 9951809


'''

import re

class alignment_obj:
    
    def __init__(self):
        print "loading and initialize the object ~~\n"
        self.seqs = {}
        self.alignment_array = []
        self.alignment_array_up = []
        self.alignment_array_left = []
        self.alignment_arg = []
        self.alignment_arg_up = []
        self.alignment_arg_left = []
        self.align_x = ""
        self.align_y = ""
        self.lcs = ""
        self.gap_open = 0
        self.gap_extend = 0
        self.amino_acid = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*');
        self.score = {}
        
    def load_seq(self,filename):    
        # finding longest common sequences
        read = open(filename,'r')
        count = 0
        for line in read:
            line = line.lstrip(' ')
            line = line.strip('\n')
            line = line.rstrip(' ')
            self.seqs[count] = line
            count = count+1

    def load_matrix(self,filename_m):
        read = open(filename_m,'r')
        count = 0;
        
        #initialize the scoring schema
        
        for i in range(0,len(self.amino_acid)):
            self.score[self.amino_acid[i]]={}
            for j in range(0,len(self.amino_acid)):
                self.score[self.amino_acid[i]][self.amino_acid[j]]=0    
            
        for line in read:
            if (line[0]=="#"):
                next
            else:
                line.lstrip(' ')
                line.strip('\n')
                line.rstrip(' ')
                a = re.split('\s+',line)
                if (count != 0):
                    temp_aa = a[0]
                    #print len(a)
                    for i in range(1,len(a)-1):
                        self.score[temp_aa][self.amino_acid[i-1]]=a[i]
                        #print "%s %s %s %s" % (i, self.score[temp_aa][self.amino_acid[i-1]],temp_aa,self.amino_acid[i-1])
                count = count+1
                
    def get_lcs(self,index_x,index_y):
        
        print "start to get the longest common sequence ~~ "
        self.align_x = ""
        self.align_y = ""
        for i in range(0,len(self.seqs[index_x])+1):
            self.alignment_array.append([])
            self.alignment_arg.append([])
            for j in range(0,len(self.seqs[index_y])+1):
                self.alignment_array[i].append('NA')
                self.alignment_arg[i].append('NA')
        
        #initialize the condition of longest common sequence
        
        for i in range(0,len(self.seqs[index_x])+1):
            self.alignment_array[i][0] = 0
        
        for j in range(0,len(self.seqs[index_y])+1):
            self.alignment_array[0][j] = 0
            
        for i in range(1,len(self.seqs[index_x])+1):
            for j in range(1,len(self.seqs[index_y])+1):
                print self.alignment_array[i][j]
            #determine who is the maximal value
                temp_max = -999999
                temp_dir = ""
                sigma = 0
                if (self.seqs[index_x][i-1] == self.seqs[index_y][j-1]):
                    #print " %s %s " % (self.seqs[index_x][i-1],self.seqs[index_y][j-1])
                    sigma = 1
                    if (int(self.alignment_array[i-1][j-1])+sigma > temp_max):
                        temp_max = int(self.alignment_array[i-1][j-1])+sigma
                        temp_dir = "ul"
                else:
                    sigma = 0
                    
                
                if (int(self.alignment_array[i-1][j]) > temp_max):
                    temp_max = int(self.alignment_array[i-1][j])
                    temp_dir = "up"
                if (int(self.alignment_array[i][j-1]) > temp_max):
                    temp_max = int(self.alignment_array[i][j-1])
                    temp_dir = "left"  
                   
                
                self.alignment_array[i][j] = temp_max
                self.alignment_arg[i][j] = temp_dir
        
        #check the alignment result
        temp = ""
        for i in range(0,len(self.seqs[index_x])):
            #print "\n"
            temp = ""
            for j in range(0,len(self.seqs[index_y])):
                temp += " %s " % self.alignment_array[i][j]
            print temp
        temp = ""                 
        for i in range(0,len(self.seqs[index_x])):
            #print "\n"
            temp = ""
            for j in range(0,len(self.seqs[index_y])):        
                temp += " %s " % (self.alignment_arg[i][j])
            print temp
    
    def get_align_aff(self,index_x,index_y,gap_op,gap_exd,semi_opt):
        
        print "start to get the alignment ~~ "
        self.align_x = ""
        self.align_y = ""
        self.gap_extend = gap_exd
        self.gap_open = gap_op
        
        #initialize the matrix
        
        for i in range(0,len(self.seqs[index_x])+1):
            self.alignment_array.append([])
            self.alignment_array_left.append([])
            self.alignment_array_up.append([])
            self.alignment_arg.append([])
            self.alignment_arg_up.append([])
            self.alignment_arg_left.append([])
            for j in range(0,len(self.seqs[index_y])+1):
                self.alignment_array[i].append('NA')
                self.alignment_array_up[i].append('NA')
                self.alignment_array_left[i].append('NA')
                self.alignment_arg[i].append('\\')
                self.alignment_arg_up[i].append('|')
                self.alignment_arg_left[i].append('-')
                
        
        #initialize the original position
        
        self.alignment_array[0][0] = 0
        #self.alignment_array_up[0][0] = 0
        #self.alignment_array_left[0][0] = 0
        
        #initialize the boundary condition of alignment (i,0) (0,j)
        
        for i in range(1,len(self.seqs[index_x])+1):
            #self.alignment_array[i][0] = -999999
            if (semi_opt == 0):
                self.alignment_array[i][0] = -999999
                self.alignment_array_up[i][0] = -self.gap_open-(i-1)*self.gap_extend
                self.alignment_array_left[i][0] = -999999
            else:
                self.alignment_array[i][0] = -self.gap_open-(i)*self.gap_extend
                #self.alignment_array_left[i][0] = -self.gap_open-(i)*self.gap_extend
                self.alignment_array_up[i][0] = -self.gap_open-(i)*self.gap_extend
            
            
            #self.alignment_array_left[i][0] = 0
        
        for j in range(1,len(self.seqs[index_y])+1):
            #self.alignment_array[0][j] = -999999
            if (semi_opt == 0):
                self.alignment_array[0][j] = -999999
                self.alignment_array_left[0][j] = -self.gap_open-(i-1)*self.gap_extend
                self.alignment_array_up[0][j] = -999999
            else:
                self.alignment_array[0][j] = -self.gap_open-(j)*self.gap_extend
                self.alignment_array_left[0][j] = -self.gap_open-(j)*self.gap_extend
                #self.alignment_array_up[0][j] = -self.gap_open-(j)*self.gap_extend
            
        
        
        self.dump_matrix(self.alignment_array_up)    
        # start the dynamic programming part
            
        for i in range(1,len(self.seqs[index_x])+1):
            for j in range(1,len(self.seqs[index_y])+1):
                #print self.alignment_array[i][j]
            #determine who is the maximal value
                print "%s %s %s %s" % (i,j,self.seqs[index_x][i-1],self.seqs[index_y][j-1])
                temp_max = -999999
                temp_max1 = -999999
                temp_max2 = -999999
                temp_dir = ""
                matched_score = float(self.score[self.seqs[index_x][i-1]][self.seqs[index_y][j-1]])
                
                # I_up(i,j) = max(self.alignment[i-1)
                print "%d %d " % (i,j)
                if (float(self.alignment_array[i-1][j])-(self.gap_open)-self.gap_extend > temp_max1):
                    self.alignment_array_left[i][j] = float(self.alignment_array[i-1][j])-(self.gap_open)-self.gap_extend
                    temp_max1 = self.alignment_array_left[i][j]
                    self.alignment_arg_left[i][j] = "M"
                
                if (float(self.alignment_array_left[i-1][j])-self.gap_extend > temp_max1):
                    self.alignment_array_left[i][j] = float(self.alignment_array_left[i-1][j])-self.gap_extend
                    self.alignment_arg_left[i][j] = "-"   
                
                    
                if (float(self.alignment_array[i][j-1])-(self.gap_open)-self.gap_extend > temp_max2):
                    self.alignment_array_up[i][j] = float(self.alignment_array[i][j-1])-(self.gap_open)-self.gap_extend
                    temp_max2 = self.alignment_array_up[i][j]
                    self.alignment_arg_up[i][j] = "M"
                
                if (float(self.alignment_array_up[i][j-1]) - self.gap_extend > temp_max2):
                    self.alignment_array_up[i][j] = float(self.alignment_array_up[i][j-1])-self.gap_extend
                    self.alignment_arg_up[i][j] = "|" 
                
                # take maximal value of max(M(i-1,j-1)+I_up(i-1,j-1)+I_left(i-1,j-1)
                
                if (float(self.alignment_array[i-1][j-1]) +matched_score > temp_max):
                    temp_max = float(self.alignment_array[i-1][j-1]) +matched_score
                    temp_dir = "\\"
                if (float(self.alignment_array_up[i][j]) > temp_max):
                    temp_max = float(self.alignment_array_up[i][j])
                    temp_dir = "U"
                if (float(self.alignment_array_left[i][j]) > temp_max):
                    temp_max = float(self.alignment_array_left[i][j])
                    temp_dir = "L"                     
                
                
                                     
                self.alignment_array[i][j] = temp_max
                self.alignment_arg[i][j] = temp_dir
        
        #check the alignment result
        temp = ""
        for i in range(0,len(self.seqs[index_x])):
            #print "\n"
            temp = ""
            for j in range(0,len(self.seqs[index_y])):
                temp += " %s " % self.alignment_array[i][j]
            #print temp
        temp = ""                 
        for i in range(0,len(self.seqs[index_x])):
            #print "\n"
            temp = ""
            for j in range(0,len(self.seqs[index_y])):        
                temp += " %s " % (self.alignment_arg[i][j])
            #print temp
        self.dump_matrix(self.alignment_array_up)
        self.dump_matrix(self.alignment_array_left)
        self.dump_matrix(self.alignment_array)
        self.dump_matrix(self.alignment_arg)
        self.dump_matrix(self.alignment_arg_up)
        self.dump_matrix(self.alignment_arg_left)
        
        print "Score: %s " % temp_max
            
    def print_lcs(self,index_x,index_y,i,j):
        
        if (i == 0 or j == 0):
            print "now printing ..."
            if (i == 0):
                #print "gag"
                for x in range(j,0,-1):
                    print "hah %s " % self.seqs[index_y][x-1]
                    self.align_x = "-" + self.align_x
                    self.align_y = self.seqs[index_y][x-1] + self.align_y.strip("\n")
            if(j == 0):
                #print "gag %d" % i
                for x in range(i,0,-1):
                    print "hah %s " % self.seqs[index_x][x-1]
                    self.align_x = self.seqs[index_x][x-1]+ self.align_x.strip("\n")
                    self.align_y = "-" + self.align_y    
            return True
        if (self.alignment_arg[i][j] == "ul"):            
            self.align_x = self.seqs[index_x][i-1] + self.align_x
            self.align_y = self.seqs[index_y][j-1] + self.align_y
            self.lcs = self.seqs[index_x][i-1]+self.lcs
            self.print_lcs(index_x,index_y,i-1,j-1)
        else:
            if (self.alignment_arg[i][j] == "up"):
                self.align_x = self.seqs[index_x][i-1]+self.align_x
                self.align_y = "-"+self.align_y
                self.print_lcs(index_x,index_y,i-1,j)
            else:
                self.align_x = "-"+self.align_x
                self.align_y = self.seqs[index_y][j-1]+self.align_y
                self.print_lcs(index_x,index_y,i,j-1)
 
    def print_aff(self,index_x,index_y,i,j,M_opt):
        if (i < 1 or j < 1): return True
        #print "%d %d %s" % (i,j,self.alignment_arg[i][j])
        if (M_opt == "M"):
            if (self.alignment_arg[i][j] == "L"):
                self.print_aff(index_x,index_y,i,j,"L")
            elif (self.alignment_arg[i][j] == "U"):
                self.print_aff(index_x,index_y,i,j,"U")
            else :
                self.align_x = self.seqs[index_x][i-1] + self.align_x
                self.align_y = self.seqs[index_y][j-1] + self.align_y
                self.print_aff(index_x,index_y,i-1,j-1,"M")
        elif (M_opt == "L"):
            if (self.alignment_arg_left[i][j] == "M"):
                self.align_x = self.seqs[index_x][i-1] + self.align_x
                self.align_y = "-" + self.align_y
                self.print_aff(index_x,index_y,i-1,j,"M")
            else:
                self.align_x = self.seqs[index_x][i-1] + self.align_x
                self.align_y = "-" + self.align_y
                self.print_aff(index_x,index_y,i-1,j,"L")
        else:
            if (self.alignment_arg_up[i][j] == "M"):
                self.align_x = "-" + self.align_x
                self.align_y = self.seqs[index_y][j-1] + self.align_y
                self.print_aff(index_x,index_y,i,j-1,"M")
            else:
                self.align_x = "-" + self.align_x
                self.align_y = self.seqs[index_y][j-1] + self.align_y
                self.print_aff(index_x,index_y,i,j-1,"U")
 
    def dump_matrix(self,a):
        text = ""
        for i in range(0,len(a)):
            text += ""
            for j in range(0,len(a[i])):
                text += "%2s " % a[i][j]
            text +="\n"
        print text
            
    def reset(self):
        print "reset all value" 
        self.alignment_array = []
        self.alignment_arg = []
        self.align_x = ""
        self.align_y = ""
        self.lcs = "" 
        
    def get_length(self,index):
        return len(self.seqs[index])
    
    def get_align_re(self):
        print "%s\n%s\n" % (self.align_x,self.align_y) 
        return (self.align_x,self.align_y)
    
    def get_lcs_re(self):
        print "Longest common sequence : %s \n" % self.lcs
        return self.lcs
        
    def reset_seq(self):
        self.seqs = {}
    
    def get_seq(self):
        return self.seqs
    
    def get_alignment_array(self):
        return self.alignment_array
    
    def get_alignment_arg(self):
        return self.alignment_arg
    
    def check_score(self):
        score_tmp = 0
        for i in range(0,len(self.align_x)):
            if (self.align_x[i] == "-" or self.align_y[i] == "-"):
                if (self.align_x[i-1] != "-" and self.align_y[i-1] != "-"):
                    score_tmp += -10-0.5
                else:
                    score_tmp += -0.5
            else:
                score_tmp += int(self.score[self.align_x[i]][self.align_y[i]]) 
        print "checked score : %s" %  score_tmp
            