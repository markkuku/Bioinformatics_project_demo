'''
Created on 2011/11/12

@author: markku

# console mode of motif finding algorithm

'''

from time import clock, time
from motif_find_obj import *
from Tkinter import *
import math
import random
import tkMessageBox 
import tkFileDialog

start_time = time()
mf = motif_finding()

def run():
    print "Running the program"
    begin_time = time()
    mf.reset_all()
    length = entry_length.get()
    Status_label.config(text="STATUS: RUNNING ...")
    print "again"
    print length
    if (length is None):
        tkMessageBox.showinfo("Message", "You shoule input the number.")
        Status_label.config(text="STATUS: ERROR !! ...")
        return False     
    else:
        try :
            int(length)
            if (length > 0):
                length=int(length)
                mf.set_motif_len(length)             
        except ValueError:
                tkMessageBox.showinfo("Message", "You shoule input the number.")
                Status_label.config(text="STATUS: ERROR !! ...")
                return False
    print "again"
    print length
    mf.set_motif_len(length) 
    best = mf.motif_finding_bfw()
    end_time = time()
    Status_label.config(text="STATUS : FINISHED !! ...")
    Location_label.config(text="Location : %s " % mf.get_locat())
    Motif_label.config(text="Best Motif: %s " % best)
    Dis_label.config(text="Best Distance: %s " % mf.get_min_dis())
    time_label.config(text="Time elapsed: %0.3f seconds" % (end_time-begin_time))
    print best
        
    #label the sequence
    
    final_seq = mf.get_seq()
    final_locat = mf.get_locat()    
    final_locat = final_locat.rstrip(']')
    final_locat = final_locat.lstrip('[')
    locat = final_locat.split(',')
    #text.delete(1.0, END)
    load_to_text(final_seq)
    for i in range(0,len(final_seq)):
        header_text = "sequence %d : " % i
        #text.insert(END,header_text)
        shift = len(header_text)
        start_l = int(locat[i])
        #text.insert(END,final_seq[i][0:start_l])
        temp_name = final_seq[i][start_l:(start_l+length)]
        text.tag_add(temp_name, "%d.%d" %((i+1),start_l+shift),"%d.%d" %((i+1),(start_l+length+shift)))
        text.tag_config(temp_name,background="yellow",foreground="red",underline=1)
        #text.insert(END,final_seq[i][(start_l+length):len(final_seq[i])])
        #text.insert(END,"\n")
        print "%d %s %d %d " % (shift,temp_name,start_l,length)

def openFileDialog(*arg):
    print "haha opening program ... !"
    mf.reset_seq()    
    filename = openFile.show()
    print filename
    mf.load_seq(filename)
    seq = mf.get_seq()
    load_to_text(seq)

def load_to_text(in_seq):
    text.delete(1.0,END)
    for i in range(0,len(in_seq)):
        text.insert(END,"sequence %d : %s \n" % (i,in_seq[i]))
    

def show_version(*arg):
    print "Version 1.00"
    tkMessageBox.showinfo("Message", "Version 1.00") 

def show_help(*arg):
    tkMessageBox.showinfo("Message", "Version 1.00")

win = Tk()
win.title('Motif Finding Program')
win.geometry("740x400")

# for the frame one

f = Frame(win)
f1 = Frame(win)

# Main menu

mainmenu = Menu(win)
filemenu = Menu(mainmenu)
filemenu.add_command(label="Open Sequence file",command=openFileDialog,accelerator="Ctrl-O")
filemenu.add_separator()
filemenu.add_command(label="Exit",command=win.quit)
win.bind("<Control-o>", openFileDialog)

#help menu

helpmenu = Menu(mainmenu)
helpmenu.add_command(label="version", command=show_version)
helpmenu.add_command(label="help", command=show_help)

mainmenu.add_cascade(label="File", menu=filemenu)
mainmenu.add_cascade(label="Help", menu=helpmenu)

# Frame 1: 

Status_label = Label(f)
Status_label.grid(row=0,column=3,pady=10,padx=10)
Status_label.config(text = "Status : N/A")
Location_label = Label(f)
Location_label.grid(row=0,column=4,pady=10,padx=10)
Location_label.config(text = "Location : N/A")
Motif_label = Label(f)
Motif_label.grid(row=0,column=5,pady=10,padx=10)
Motif_label.config(text = "Best Motif : N/A")
Dis_label = Label(f)
Dis_label.grid(row=1,column=4,pady=10,padx=10)
Dis_label.config(text = "Best Distance : N/A")
time_label = Label(f)
time_label.grid(row=1,column=5,pady=10,padx=10)
time_label.config(text = "Elapse time : N/A")
enter_n = Button(f, text = "Run", command = run)
enter_n.grid(row=0,column=2,pady=10,padx=10)
show_length = Label(f)
show_length.config(text = "Motif Length:")
show_length.grid(row=0,column=0,pady=10,padx=10)
entry_length = Entry(f,width = 3)
entry_length.insert(0,"3")
entry_length.grid(row=0,column=1,pady=10,padx=10)
openFile = tkFileDialog.Open(win)

# Frame 2:

s = Scrollbar(f1)
text = Text(f1)
text.config(width = 100, height = 20)
s.pack(side = RIGHT,fill=Y)
text.pack(side = LEFT,fill=Y)
s.config(command=text.yview)

text.config(yscrollcommand=s.set)

f1.pack(side = BOTTOM, fill=BOTH)
f.pack(side = TOP, fill = BOTH)
win.config(menu=mainmenu)
win.mainloop()
