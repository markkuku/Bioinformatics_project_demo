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
        self.alignment_array_M = []
        self.alignment_array_up = []
        self.alignment_array_left = []
        self.alignment_arg = []
        self.align_x = ""
        self.align_y = ""
        self.lcs = ""
        self.gap_open = 0
        self.gap_extend = 0
        self.amino_acid = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*');
        self.score = {}
        self.x = 0
        self.y = 1
        
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
    
    def get_align_aff(self,index_x,index_y,gap_op,gap_exd,semi_opt,time_delay,win,textfield,print_opt):
        
        print "start to get the alignment ~~ %s" % time_delay
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
            self.alignment_array_M.append([])
            for j in range(0,len(self.seqs[index_y])+1):
                self.alignment_array[i].append('NA')
                self.alignment_array_M[i].append('NA')
                self.alignment_array_up[i].append('NA')
                self.alignment_array_left[i].append('NA')
                self.alignment_arg[i].append('NA')
        
        #initialize the original position
        self.alignment_array[0][0] = 0
        self.alignment_array_M[0][0] = 0
        self.alignment_array_up[0][0] = 0
        self.alignment_array_left[0][0] = 0
        
        #initialize the boundary condition of alignment (i,0) (0,j)
        
        for i in range(1,len(self.seqs[index_x])+1):
            #self.alignment_array[i][0] = -999999
            if (semi_opt == 0):
                self.alignment_array[i][0] = -999999
                self.alignment_array_M[i][0] = -999999
                self.alignment_array_up[i][0] =  -self.gap_open-(i-1)*self.gap_extend
                self.alignment_array_left[i][0] = -999999
            else:
                self.alignment_array_M[i][0] = -self.gap_open-(i-1)*self.gap_extend
                self.alignment_array[i][0] = -self.gap_open-(i-1)*self.gap_extend
                #self.alignment_array_left[i][0] = -self.gap_open-(i-1)*self.gap_extend
                self.alignment_array_up[i][0] = -self.gap_open-(i-1)*self.gap_extend
            
            
            #self.alignment_array_left[i][0] = 0
        
        for j in range(1,len(self.seqs[index_y])+1):
            #self.alignment_array[0][j] = -999999
            if (semi_opt == 0):
                self.alignment_array[0][j] = -999999
                self.alignment_array_M[0][j] = -999999
                self.alignment_array_left[0][j] = -self.gap_open-(j-1)*self.gap_extend
                self.alignment_array_up[0][j] =  -999999
            else:
                self.alignment_array_M[0][j] = -self.gap_open-(j-1)*self.gap_extend
                self.alignment_array[0][j] = -self.gap_open-(j-1)*self.gap_extend
                self.alignment_array_left[0][j] = -self.gap_open-(j-1)*self.gap_extend
                #self.alignment_array_up[0][j] = -self.gap_open-(j-1)*self.gap_extend
            
        
        
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
                self.x = index_x
                self.y = index_y
                matched_score = float(self.score[self.seqs[index_x][i-1]][self.seqs[index_y][j-1]])
                #M(i,j) = V(i-1,j-1) + m/mm
                self.alignment_array_M[i][j] = self.alignment_array[i-1][j-1]+matched_score
                self.print_proc(i, j, 4, 4, 2, time_delay, win, textfield,print_opt)
                if (float(self.alignment_array_M[i-1][j])-(self.gap_open) > temp_max1):
                    self.alignment_array_left[i][j] = float(self.alignment_array_M[i-1][j])-(self.gap_open)
                    temp_max1 = self.alignment_array_left[i][j]
                    self.print_proc(i, j, 4, 4, 2, time_delay, win, textfield,print_opt)
                if (float(self.alignment_array_left[i-1][j])-self.gap_extend > temp_max1):
                    self.alignment_array_left[i][j] = float(self.alignment_array_left[i-1][j])-self.gap_extend   
                    self.print_proc(i, j, 4, 4, 2, time_delay, win, textfield,print_opt)
                
                self.print_proc(i, j, 4, 4, 3, time_delay, win, textfield,print_opt)    
                if (float(self.alignment_array_M[i][j-1])-(self.gap_open) > temp_max2):
                    self.alignment_array_up[i][j] = float(self.alignment_array_M[i][j-1])-(self.gap_open)
                    temp_max2 = self.alignment_array_up[i][j]
                    self.print_proc(i, j, 4, 4, 3, time_delay, win, textfield,print_opt)
                
                if (float(self.alignment_array_up[i][j-1]) - self.gap_extend > temp_max2):
                    self.alignment_array_up[i][j] = float(self.alignment_array_up[i][j-1])-self.gap_extend
                    self.print_proc(i, j, 4, 4, 3, time_delay, win, textfield,print_opt)   
                self.print_proc(i, j, 4, 4, 3, time_delay, win, textfield,print_opt)     
                
                # take maximal value of max(M(i-1,j-1)+I_up(i-1,j-1)+I_left(i-1,j-1)
                self.print_proc(i, j, 4, 4, 1, time_delay, win, textfield,print_opt)
                if (float(self.alignment_array_M[i][j]) > temp_max):
                    temp_max = float(self.alignment_array_M[i][j])
                    temp_dir = "ul"
                    self.print_proc(i, j, 4, 4, 1, time_delay, win, textfield,print_opt)
                if (float(self.alignment_array_up[i][j]) > temp_max):
                    temp_max = float(self.alignment_array_up[i][j])
                    temp_dir = "up"
                    self.print_proc(i, j, 4, 4, 1, time_delay, win, textfield,print_opt)
                if (float(self.alignment_array_left[i][j]) > temp_max):
                    temp_max = float(self.alignment_array_left[i][j])
                    temp_dir = "left"
                    self.print_proc(i, j, 4, 4, 1, time_delay, win, textfield,print_opt)                     
                self.alignment_array[i][j] = temp_max
                self.alignment_arg[i][j] = temp_dir
                self.print_proc(i, j, 4, 4, 1, time_delay, win, textfield,print_opt)
                # I_up(i,j) = max(self.alignment[i-1)
                               
               
        

        self.dump_matrix(self.alignment_array_up)
        self.dump_matrix(self.alignment_array_left)
        self.dump_matrix(self.alignment_array)
        self.dump_matrix(self.alignment_arg)
        
        print "Score: %s " % temp_max
            
    def print_lcs(self,index_x,index_y,i,j,win,text,time_delay,open):
        self.print_bk(i,j,4,4,win,text,time_delay,open)
        if (i == 0 or j == 0):
            print "now printing ..."
            if (i == 0):
                #print "gag"
                for x in range(j,0,-1):
                    #print "hah %s " % self.seqs[index_y][x-1]
                    self.align_x = "-" + self.align_x
                    self.align_y = self.seqs[index_y][x-1] + self.align_y.strip("\n")
            if(j == 0):
                #print "gag %d" % i
                for x in range(i,0,-1):
                    #print "hah %s " % self.seqs[index_x][x-1]
                    self.align_x = self.seqs[index_x][x-1]+ self.align_x.strip("\n")
                    self.align_y = "-" + self.align_y    
            return True
        if (self.alignment_arg[i][j] == "ul"):            
            self.align_x = self.seqs[index_x][i-1] + self.align_x
            self.align_y = self.seqs[index_y][j-1] + self.align_y
            self.lcs = self.seqs[index_x][i-1]+self.lcs
            self.print_lcs(index_x,index_y,i-1,j-1,win,text,time_delay,open)
        else:
            if (self.alignment_arg[i][j] == "up"):
                self.align_x = self.seqs[index_x][i-1]+self.align_x
                self.align_y = "-"+self.align_y
                self.print_lcs(index_x,index_y,i-1,j,win,text,time_delay,open)
            else:
                self.align_x = "-"+self.align_x
                self.align_y = self.seqs[index_y][j-1]+self.align_y
                self.print_lcs(index_x,index_y,i,j-1,win,text,time_delay,open)
    
    def print_proc(self,index_x,index_y,x,y,opt,time_delay,win,textfield,open):
        if (open == 0): return True
        textfield.delete(1.0, "end")
        opt_des = ['','Main recursion: \n S_main[i,j] = Max(S_Main[i-1,j-1]+BLO(Si,Sj), S_left[i-1,j-1]+BLO(Si,Sj), S_up[i-1,j-1]+BLO(Si,Sj))'
                   ,'Left recursion: \n S_Left[i,j] = Max(S_Main[i-1,j]-d-e, S_Left[i-1,j]-e)'
                   ,'up recursion: \n S_Up[i,j] = Max( S_main[i,j-1]-d-e, S_Up[i,j-1]-e)']
        m_desc = ['','S_main[i,j]','S_Left[i,j]','S_Up[i,j]']
        textfield.insert("end",opt_des[opt]+"\n\n") 
        range_print_x = int(x/2)
        range_print_y = int(y/2)
        upper_b_x = len(self.seqs[0])
        upper_b_y = len(self.seqs[1])
        
        text_tmp = ""
        if ((index_x - range_print_x) < 0):
            x_lower = 0
            x_upper = x 
        elif ((index_x + range_print_x)> upper_b_x):
            x_lower = upper_b_x - x + 1 
            x_upper = upper_b_x + 1
        else:
            x_lower = index_x - range_print_x
            x_upper = index_x + range_print_x
            
        if ((index_y - range_print_y) < 0):
            y_lower = 0
            y_upper = y
        elif ((index_y + range_print_y) > upper_b_y):
            y_lower = upper_b_y - y + 1
            y_upper = upper_b_y + 1
        else:
            y_lower = index_y - range_print_y
            y_upper = index_y + range_print_y 
        print "%d %d %d %d %d %d \n" % (index_x, index_y,x_lower ,x_upper,y_lower, y_upper ) 
        text_tmp += "    "
        
        print_get_max = [0,0,0]
        matched_score = float(self.score[self.seqs[self.x][index_x-1]][self.seqs[self.y][index_y-1]])
        if (opt == 1):
            print_get_max[0] = self.alignment_array[index_x-1][index_y-1]+matched_score
            print_get_max[1] = self.alignment_array_up[index_x-1][index_y-1]+matched_score
            print_get_max[2] = self.alignment_array_left[index_x-1][index_y-1]+matched_score
            text_tmp += "Current: Max ( %s, %s, %s ) Direction: (upper left, upper, left)\n" % (print_get_max[0],print_get_max[1],print_get_max[2]) 
        elif (opt == 2):
            print_get_max[0] = self.alignment_array[index_x-1][index_y]-self.gap_extend-self.gap_open
            print_get_max[1] = self.alignment_array_left[index_x-1][index_y]-self.gap_extend
            text_tmp += "Current: Max ( %s, %s ) \n" % (print_get_max[0],print_get_max[1]) 
        elif (opt == 3):
            print_get_max[0] = self.alignment_array[index_x][index_y-1]-self.gap_extend-self.gap_open
            print_get_max[1] = self.alignment_array_up[index_x][index_y-1]-self.gap_extend
            text_tmp += "Current: Max ( %s, %s ) \n" % (print_get_max[0],print_get_max[1])
        text_tmp += "\n   "
        if (opt == 3):
            text_tmp += " %28s       %28s" % (m_desc[1],m_desc[3])
        elif (opt == 2):
            text_tmp += " %28s       %28s" % (m_desc[1],m_desc[2])
        else:
            text_tmp += " %28s       %28s       %28s" % (m_desc[1],m_desc[2],m_desc[3])
        text_tmp += "       "   
        text_tmp += "\n\n   "
        for k in range(0,3):
            for i in range(x_lower,x_upper):
                text_tmp += " %7s" % i
            text_tmp += "     "
        text_tmp += "\n"
        
        for j in range(y_lower,y_upper):                        
            text_tmp += " %d " % j
            for i in range(x_lower,x_upper):
                
                if (i == index_x - 1 and j == index_y -1 and opt == 1):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x -1 and j == index_y and opt == 2):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x and j == index_y-1 and opt == 3):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x and j == index_y and opt == 1):
                    temp_array = "["+self.ch_to_inf(self.alignment_array[i][j])+"]"
                    text_tmp += " %7s" % temp_array
                else:
                    text_tmp += " %7s" % self.ch_to_inf(self.alignment_array[i][j])
            text_tmp += "     "
            for i in range(x_lower,x_upper):
                if (opt == 3):
                    continue;
                if (i == index_x - 1 and j == index_y -1 and opt == 1):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array_left[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x -1 and j == index_y and opt == 2):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array_left[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x and j == index_y and opt == 2):
                    temp_array = "["+self.ch_to_inf(self.alignment_array_left[i][j])+"]"
                    text_tmp += " %7s" % temp_array
                else:
                    text_tmp += " %7s" % self.ch_to_inf(self.alignment_array_left[i][j])
            if (opt != 3):
                text_tmp += "     "   
            for i in range(x_lower,x_upper):
                if (opt == 2):
                    continue;
                if (i == index_x - 1 and j == index_y -1 and opt == 1):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array_up[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x  and j == index_y -1 and opt == 3):
                    temp_array = "<"+self.ch_to_inf(self.alignment_array_up[i][j])+">"
                    text_tmp += " %7s" % temp_array
                elif (i == index_x and j == index_y and opt == 3):
                    temp_array = "["+self.ch_to_inf(self.alignment_array_up[i][j])+"]"
                    text_tmp += " %7s" % temp_array
                else:
                    text_tmp += " %7s" % self.ch_to_inf(self.alignment_array_up[i][j])                     
            text_tmp += "\n"
        text_tmp += "\n"
        for j in range(y_lower,y_upper):                        
            text_tmp += " %d " % j
            for i in range(x_lower,x_upper):               
                if (i == index_x and j == index_y):
                    temp_array = "["+self.ch_to_inf(self.alignment_arg[i][j])+"]"
                    text_tmp += " %7s" % temp_array
                else:
                    text_tmp += " %7s" % self.ch_to_inf(self.alignment_arg[i][j])    
            text_tmp += "\n"
        textfield.insert("end",text_tmp)
        #time.sleep(0.25)
        #print time_delay
        win.after(time_delay)
        textfield.update_idletasks()

    def print_bk(self,index_x,index_y,x,y,win,textfield,time_delay,open):
        if (open == 0):
            return True
        textfield.delete(1.0, "end")
        textfield.insert("end","Now Backtracking process ............ \n\n")
        range_print_x = int(x/2)
        range_print_y = int(y/2)
        upper_b_x = len(self.seqs[0])
        upper_b_y = len(self.seqs[1])
        text_tmp = ""
        if ((index_x - range_print_x) < 0):
            x_lower = 0
            x_upper = x 
        elif ((index_x + range_print_x)> upper_b_x):
            x_lower = upper_b_x - x + 1 
            x_upper = upper_b_x + 1
        else:
            x_lower = index_x - range_print_x
            x_upper = index_x + range_print_x 
            
        if ((index_y - range_print_y) < 0):
            y_lower = 0
            y_upper = y
        elif ((index_y + range_print_y) > upper_b_y):
            y_lower = upper_b_y - y + 1
            y_upper = upper_b_y + 1
        else:
            y_lower = index_y - range_print_y
            y_upper = index_y + range_print_y  
        text_tmp += "    "
       
        for j in range(y_lower,y_upper):                        
            text_tmp += " %d " % j
            for i in range(x_lower,x_upper):               
                if (i == index_x and j == index_y):
                    temp_array = "["+self.ch_to_inf(self.alignment_arg[i][j])+"]"
                    text_tmp += " %7s" % temp_array
                else:
                    text_tmp += " %7s" % self.ch_to_inf(self.alignment_arg[i][j])    
            text_tmp += "\n"
        text_tmp += "\n"
        text_tmp += "Current alignment: \n %s\n %s\n \n\n" % (self.align_x,self.align_y)
        textfield.insert("end",text_tmp)
        win.after(time_delay)
        textfield.update_idletasks()
    
    def dump_matrix(self,a):
        text = ""
        for i in range(0,len(a)):
            text += ""
            for j in range(0,len(a[i])):
                text += "%2s " % a[i][j]
            text +="\n"
        print text
    
    def ch_to_inf(self,num):
        if (num == -999999):
            return "-inf"
        else:
            return str(num)
            

        
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
    
    def get_alignment_array_up(self):
        return self.alignment_array_up
    
    def get_alignment_array_left(self):
        return self.alignment_array_left
    
    def get_alignment_arg(self):
        return self.alignment_arg
    
    def get_score(self):
        return self.alignment_array[len(self.seqs[self.x])][len(self.seqs[self.y])]
    
    def reset(self):
        self.alignment_array = []
        self.alignment_array_M = []
        self.alignment_array_up = []
        self.alignment_array_left = []
        self.alignment_arg = []
        self.align_x = ""
        self.align_y = ""
        self.lcs = ""
        self.score = {}
        