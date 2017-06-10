Homework Project IV: Affine gap alignment with BLOSUM 62 (printing procedure is implemented) 
Name: Shih-Yen Ku
TIGP ID: 9951809

File description:
1. alignment_obj.py: contains the class and functions
2. alignment_console.py: contains the console mode of alignment program (obsolete I think)
3. align_prj4_gui.py: contains the GUI mode of alignment program (include the printing procedure)

Input Example: (two sequence)
TGCATA
ATCTGAT

Usage:
1. Console mode: python alignment_console.py (obsolete)
2. GUI mode: double click the motif_prj3_gui.py

Output:
Alignment result
The first matrix is score
The second matrix is storing the directions. ul: upper left; left: left; up: up.

Example.

Score: -3.0 

NGPS
QNQL

  0 -999999 -999999 -999999 -999999 
-999999 0.0 -5.5 -11.5 -14.5 
-999999 -12.0 0.0 -7.5 -15.0 
-999999 -11.5 -12.5 -1.0 -10.5 
-999999 -11.0 -10.0 -10.5 -3.0 

  0 -999999 -999999 -999999 -999999 
-10.0 -10.5 -10.5 -11.0 -11.5 
-10.5 -11.0 -11.5 -10.5 -11.0 
-11.0 -11.5 -12.0 -12.5 -11.5 
-11.5 -12.0 -12.5 -13.0 -13.5 