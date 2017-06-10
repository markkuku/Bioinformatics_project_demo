'''
Created on 2011/12/1

@author: markkuku

ID: 9951809

'''
from alignment_obj import *

align = alignment_obj()
filename = raw_input("input your file name:\n")
gap_open = raw_input("input your gap opening penalty:\n")
gap_ext = raw_input("input your gap extension penalty:\n")
semi_opt = raw_input("Continue align with gap in the starting point or not ? (1: Yes 0: No) \n")
#filename = "test1.txt"
#gap_open = 10
#gap_ext = 0.5
#semi_opt = 0
#filename = "test2.txt"
align.load_seq(filename)
align.load_matrix("BLOSUM62.txt")
align.get_align_aff(0,1,float(gap_open),float(gap_ext),int(semi_opt))
#align.get_lcs(0,1)
align.print_aff(0, 1, align.get_length(0), align.get_length(1),"M")
align.get_align_re()
align.check_score()
#align.get_lcs_re()