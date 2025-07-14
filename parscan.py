

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
new_row = 44 #29 row;;;43 magnet coil;;;;38 gun mag;;;57 boost phase;;;65 boostersol
new_row=new_row-1
new_col = 4 #5th column, 4 is inital phase ,5 is gradient ;;;2 is the scale coil;;;; 2 is scale gun magnet

laser_idealpos=-0.0188
poslope=0.0000375116 #cm/ps, assume 0.4 eV drift
jitrange=10 #+- ;ps

#the following 3 lines define the range of phase in the scan.
nsteps =40 #from 0.1 to 10; 2.62 is 133G
inival = laser_idealpos-10*2*poslope #start from -20 ps
delval = poslope
orgfile='rr7_laserjit'
filename= orgfile+'_temp.inp'
os.system('cp '+orgfile+'.inp '+filename)
val = inival
newValueIndex = 0
while newValueIndex < nsteps:
    test5file = open(filename,'r')
    lines = test5file.readlines()
    test5file.close()

    # modifiy line 44, number 8
    new_line = string.split(lines[new_row])
   
    new_line[new_col] = str(val)
    new_line = string.join(new_line) + '\n'

    lines[new_row] = new_line
    test5file2 = open(filename, 'w')
    test5file2.writelines(lines)
    test5file2.close()
       
    os.system('parmela '+filename)
    os.system('sed -n \'49800,+0p\' TIMESTEPEMITTANCE.TBL >> tmpeng')
    os.system('echo '+str(val)+' >> tmpphase')
    #os.system('sed -n \'2500,+0p\' TIMESTEPEMITTANCE.TBL')
    #os.system('mv TIMESTEPEMITTANCE.TBL TIMESTEPEMITTANCE'+str(val)+'.TBL')
    #os.system('mv OUTPAR.TXT OUTPAR'+str(val)+'.TXT')
    #os.system('~/RTX/qplot')
    #os.system('mv emit4d.agr emit4d'+str(val)+'.agr')
    #os.system('mv emitn.agr emitn'+str(val)+'.agr')
    newValueIndex = newValueIndex+1
    print 'index: ',val
    val = val + delval
    
#newfold='lenscan'
#os.system('mkdir '+newfold)
#os.system('mv TIMES*.TBL '+newfold)
#os.system('mv emit4d*.agr '+newfold)
#os.system('mv emitn*.agr '+newfold)
#os.system('cp '+filename+' '+newfold)



os.system('paste tmpphase tmpeng > laserjit')
os.system('rm tmpeng tmpphase')
print('done')