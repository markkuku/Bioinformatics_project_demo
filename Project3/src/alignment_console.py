'''
Created on 2011/12/1

@author: markkuku

ID: 9951809

'''
from alignment_obj import *

align = alignment_obj()
filename = raw_input("input your file name:\n")
align.load_seq(filename)
align.get_lcs(0,1)
align.print_lcs(0, 1, align.get_length(0), align.get_length(1))
align.get_align_re()
align.get_lcs_re()