'''
Created on 2011/11/12

@author: markku

# console mode of motif finding algorithm

'''

from time import clock, time
from alignment_obj import *
from Tkinter import *
import math
import random
import tkMessageBox 
import tkFileDialog


#default option is brute force

option = 0
start_time = time()
align = alignment_obj()

def run():
    print "Running the program"
    print_option = var.get()
    print print_option
    time_delay = entry_delay.get()
    Status_label.config(text="STATUS: RUNNING ...")
    begin_time = time()    
    align.reset()

    if (option == 0):
        print "longest common subsequence"
        #align.load_seq("testaa.txt")
        align.get_lcs(0,1)
        align.print_lcs(0, 1, align.get_length(0), align.get_length(1))
        lcs_re = align.get_lcs_re()
        text.insert(END,"Longest Common Sequence: %s \n" % (lcs_re))
        arr_arg = align.get_alignment_arg()
        arr_org = align.get_alignment_array()
        text_tmp = ""
        (a,b) = align.get_align_re()
        text_tmp += a+"\n"+b+"\n"
        for i in range(0,len(arr_arg)):
            for j in range (0,len(arr_arg[i])):
                text_tmp += "%2s " % arr_org[i][j]
            text_tmp += "\n"       
        text_tmp += "\n"
        for i in range(0,len(arr_arg)):
            for j in range (0,len(arr_arg[i])):            
                text_tmp += "%4s " % arr_arg[i][j]    
            text_tmp += "\n"
        text.insert(END,"%s \n" % (text_tmp))
        
    else:
        print "Alignment with BLOSUM 62"

    end_time = time()
    time_label.config(text="Time elapsed: %0.3f seconds" % (end_time-begin_time))
    Status_label.config(text="STATUS : FINISHED !! ...")    


def openFileDialog(*arg):
    print "opening program ... !"
    align.reset_seq()    
    filename = openFile.show()
    print filename
    align.load_seq(filename)
    seq = align.get_seq()
    load_to_text(seq)

def load_to_text(in_seq):
    text.delete(1.0,END)
    for i in range(0,len(in_seq)):
        text.insert(END,"sequence %d : %s \n" % (i,in_seq[i]))
    

def show_version(*arg):

    tkMessageBox.showinfo("Message", "Version 1.03") 

def show_help(*arg):
    tkMessageBox.showinfo("Message", "Version 1.03")

def set_bb(*arg):
    global option
    option = 0
    opt_label.config(text = "Current: Longest Common Subsequence")
    tkMessageBox.showinfo("Message", "Longest Common Subsequence selected!") 

def set_bnb(*arg):
    global option 
    option = 1
    opt_label.config(text = "Current: Alignment with BLOSUM 62")
    tkMessageBox.showinfo("Message", " Alignment with BLOSUM 62") 
     
win = Tk()
win.title('Motif Finding Program')
win.geometry("940x480")

# for the frame one

f = Frame(win)
f1 = Frame(win)
f2 = Frame(win)

# Main menu

mainmenu = Menu(win)
filemenu = Menu(mainmenu)
filemenu.add_command(label="Open Sequence file",command=openFileDialog,accelerator="Ctrl-O")
filemenu.add_separator()
filemenu.add_command(label="Exit",command=win.quit)
win.bind("<Control-o>", openFileDialog)

#Running type menu

typemenu = Menu(mainmenu)
typemenu.add_command(label="Longest common sequence", command=set_bb)
typemenu.add_command(label="Sequence alignment with BLOSUM62", command=set_bnb)



#help menu

helpmenu = Menu(mainmenu)
helpmenu.add_command(label="version", command=show_version)
helpmenu.add_command(label="help", command=show_help)

mainmenu.add_cascade(label="File", menu=filemenu)
mainmenu.add_cascade(label="Help", menu=helpmenu)
mainmenu.add_cascade(label="Method", menu=typemenu)


# Frame 1: 

Status_label = Label(f)
Status_label.grid(row=0,column=0,pady=10,padx=10)
Status_label.config(text = "Status : N/A")
opt_label = Label(f)
opt_label.grid(row=0,column=1,pady=10,padx=10)
opt_label.config(text = "Current : Longest common sequence")
time_label = Label(f)
time_label.grid(row=0,column=2,pady=10,padx=10)
time_label.config(text = "Elapse time : N/A")
enter_n = Button(f, text = "Run", command = run)
enter_n.grid(row=0,column=3,pady=10,padx=10)
openFile = tkFileDialog.Open(win)
var = IntVar()
check_pt =  Checkbutton(f,text = "Print Procedure", variable = var, onvalue = "1", offvalue="0")
check_pt.grid(row=0,column=4)

entry_delay = Entry(f,width = 5)
entry_delay.insert(0,"100")
entry_delay.grid(row=0,column=6,sticky=E)

unit_l = Label(f)
unit_l.config(text="milisecs delay")
unit_l.grid(row=0,column=7,sticky=W)


# Frame 2:

s = Scrollbar(f1)
text = Text(f1)
text.config(width = 128, height = 26)
#text.pack(side = LEFT,fill=Y)
text.grid(row=0,column=1,pady=1,padx=1)
s.pack(side = RIGHT,fill=Y)
text.pack(side = LEFT,fill=Y)
s.config(command=text.yview)


#text1 = Label(f1)
#text1.config(text = "print tree")
#text1.grid(row=1,column=1,pady=1,padx=1)

#text1 = Text(f1)
#text1.config(width = 128, height = 14)
#text.pack(side = LEFT,fill=Y)
#text1.grid(row=1,column=1,pady=1,padx=1)



f1.pack(side = BOTTOM)
f.pack(side = TOP)
win.config(menu=mainmenu)
win.mainloop()
