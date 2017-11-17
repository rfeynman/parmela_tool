#!/usr/bin/env python

import os, sys, string, pdb

#
# Note: python's index start from 0
# On Windows, this script has to be run under Cygwin
# The Cygwin can downloaded from: http://cygwin.com/install.html
#---------------------
# 
#Author: Erdong
#Description: Energy/bunch length/...based on parmela timestep scan of driven phase (output file engout, col. 5 vs. 1).
#           

#new_row is the line number -1 of that phase variable in the ImpactT.in file.
#new_col is the location of phase variable in that line - 1. 
import numpy as np
input_start=202
input_end=215

start_row = input_start-1#74 #
num_row=input_end-start_row#350 end
new_col = 4 #5th column, 4 is inital phase ,5 is gradient ;;;2 is the scale coil;;;; 2 is scale gun magnet
row_array=[]
for n in range(num_row):
    row_array.append(start_row)
    start_row+=1
#the following 3 lines define the range of phase in the scan.
    #coil start at 24, end on 270. 24~26: cathode; 27~34: gun; 35~44 trans
    #45~270 focusing channel
nsteps= 13 #from 0.1 to 10; 2.62 is 133G
inival =50 #;;460
diff=0#;;25
delval=10
endval=inival+diff
filename='rr3.inp'
val = inival
newValueIndex = 0
while newValueIndex < nsteps:
    test5file = open(filename,'r')
    lines = test5file.readlines()
    test5file.close()
    test5file2=open(filename,'w')
    # modifiy line 44, number 8
    k=0
    for line in row_array:        
        new_line = string.split(lines[line])
        new_line[new_col] = str(val+k*diff/num_row)
        new_line = string.join(new_line) + '\n'
        lines[line] = new_line
        k+=1
	print(row_array[0])
	stline=string.split(lines[row_array[0]-2])
	print(stline[new_col])
	stline[new_col]=str(val+90)
	stline=string.join(stline)+'\n'
	lines[row_array[0]-2]=stline
	
    test5file2.writelines(lines)
    test5file2.close()
       
    os.system('parmela '+filename)
	#os.system('echo '+str(val)+' >> tmpeng')
    #os.system('sed -n $,+0p TIMESTEPEMITTANCE.TBL >> tmpeng')
    #os.system('sed -n \'8000,+0p\' TIMESTEPEMITTANCE.TBL')
    os.system('mv TIMESTEPEMITTANCE.TBL TIMESTEPEMITTANCE'+str(val)+'.TBL')
    os.system('mv OUTPAR.TXT OUTPAR'+str(val)+'.TXT')
    #os.system('~/RTX/qplot')
    #os.system('mv emit4d.agr emit4d'+str(val)+'.agr')
    #os.system('mv emitn.agr emitn'+str(val)+'.agr')
    newValueIndex = newValueIndex+1
    val = val + delval
    print 'index: ',val
newfold='buncher'
os.system('mkdir '+newfold)
os.system('mv TIMES*.TBL '+newfold)
os.system('mv OUTPAR*.TXT '+newfold)
#os.system('mv emit4d*.agr '+newfold)
#os.system('mv emitn*.agr '+newfold)
os.system('cp '+filename+' '+newfold)


#os.system('paste tmpphase tmpeng > engout')
#os.system('rm tmpeng tmpphase')
#os.system('mv engout '+newfold)
