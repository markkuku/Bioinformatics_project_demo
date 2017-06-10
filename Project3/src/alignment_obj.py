'''
Created on 2011/12/1

author: markku

ID: 9951809


'''

class alignment_obj:
    
    def __init__(self):
        print "loading and initialize the object ~~\n"
        self.seqs = {}
        self.alignment_array = []
        self.alignment_arg = []
        self.align_x = ""
        self.align_y = ""
        self.lcs = ""
        
        
    def load_seq(self,filename):    
        # finding longest common sequences
        read = open(filename,'r')
        count = 0
        for line in read:
            line.lstrip(' ')
            line.strip('\n')
            line.rstrip(' ')
            self.seqs[count] = line
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