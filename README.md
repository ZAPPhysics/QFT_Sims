# QFT_Sims

The two main files are clftanimv_2.py and qftanim.py. When run, each of these files should create a new directory in the current 
working directory where they will place a series of .png images corresponding to each time step. The only necessary inputs should
be the number of oscillators (N), the max time to run the simulation (tmax), the starting time (t), and the starting label number
(label). These inputs are found at the bottom of each program.

For large N values, the program runs noticeably slower, so the final two inputs were included so that a single simulation can be
run over multiple sessions. The recommended way to do this is as follows: in the output directory, find the most recently created
image file (e.g. img134.png). Then, change the starting time to this number divided by 10 (e.g. t=13.4) and the starting label to this
number plus 1 (e.g. label=135). This will continue the program from where it left off without overwriting any old files.

The .pdf file contained in the repository gives a description of the math used for creating the code.
