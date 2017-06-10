'''
Created on 2011/11/23

@author: markku
'''


remaining = 0

def snooze (secs):
    from time import clock, time
    import Tkinter 
    import math
    import random
    import tkMessageBox 
    import tkFileDialog
    global remaining
    root = Tkinter.Tk()
    prompt = 'hello'
    label1 = Tkinter.Label(root, text=prompt, width=len(prompt))
    label1.pack()

    remaining = secs

    def decrement_label ():
        global remaining
        text = "Snoozing %d sec(s)" % remaining
        remaining -= 1
        label1.config(text=text, width=100)
        label1.update_idletasks()

    for i in range(1, secs + 1):
        root.after(i * 1000, decrement_label )

    root.after((i+1) * 1000, lambda : root.destroy())
    root.mainloop()

snooze(30)