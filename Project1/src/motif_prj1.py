'''
Created on 2011/11/5

@author: markku
'''

from motif_find_obj import *

mf = motif_finding()
filename = 'test'
mf.load_seq(filename)
mf.set_motif_len(3)
mf.motif_finding_bfw()
