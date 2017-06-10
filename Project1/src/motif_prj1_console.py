'''
Created on 2011/11/5

@author: markku

# The console of motif finding algorithm

'''
from time import clock, time
from motif_find_obj import *


mf = motif_finding()
filename = raw_input("input your file name:\n")
length = int(raw_input("input your motif length:\n"))
mf.load_seq(filename)
mf.set_motif_len(length)
start_time = time()
best = mf.motif_finding_bfw()
end_time = time()

e_time =  end_time - start_time  
print "%s %s" % (mf.get_locat(),mf.get_min_dis())
print "%s " % best
print "it spends %0.2f seconds " % e_time