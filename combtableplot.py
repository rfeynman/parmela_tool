#!/usr/bin/env python
'''
Created on Jul 19, 2017

@author: wange
combine multiple colume into sigle TIMESTEP file. choose the coleum and phase steps and run the code
'''

import numpy as np
from numpy import genfromtxt, dtype, loadtxt
import string



fold = "test"
term = "TIMESTEPEMITTANCE"
suf = ".TBL"
begin_number=90
end_number=132
intvel=3
listrange = list(range(90, end_number+1, intvel)) # generate the number list of the TIMESTEPEMITT???
print(listrange)

def linepick(pick):
    sumdata=[]
    datanum=[]
    for intr in listrange[1:]:#iteration for open multiple TIMESTEP file
        with open("../" + fold + "/" + term + str(intr) + suf) as f:
            l = [line.split() for line in f]
        f.close()
        usefuldata = l[84:] #84 is starting the table header

        # content=genfromtxt("../"+fold+"/"+term+"90"+suf,delimiter='/t',names=True,dtype=None)
        headname = [elem for elem in usefuldata[0]][1:]
        usefuldata[0]=headname
        col_usefuldata={usefuldata[0][i]:[line[i] for line in usefuldata] for i in range(len(usefuldata[0]))} # generate the row as dictionary and rotate{header:head xx xx... xx}
        for i in range(len(headname)):
            col_usefuldata[headname[i]][0]=headname[i]+str(intr) # add number on the header
            
        for elem in pick:
            sumdata.append(col_usefuldata[headname[elem]]) #generate a table with multiple lines 
        #print(col_usefuldata[headname[0]])
        datanum.append(len(usefuldata)) #check the multiple files length, only minimum list will be applied
    minnum=min(datanum)
    with open("../" + fold + "/" + term + str(listrange[0]) + suf) as f0: #open first file
        l0=f0.readlines()
    #print(l0[84])
    f0.close()
    
    #print(l0[84].strip('\n'),str(sumdata[0][0]),"\n",'    '.join([l0[84].strip('\n'),str(sumdata[0][0])]))

    
    for k in range(84,minnum):
        for i in range(len(sumdata)):
            
            l0[k]='    '.join([l0[k].strip('\n'),str(sumdata[i][k-84])+'\n']) # append the new table to the end
            
    inslabel=l0.index('Ezref(MV/m)\n')
    
    del l0[(minnum):]

    for i in range(len(sumdata)): #insert the title and label
        l0.insert(inslabel+i+1,str(sumdata[i][0]+'\n'))
        l0.insert(inslabel+2*i+25,str(sumdata[i][0]+'\n'))
        l0.insert(inslabel+3*i+49,str(sumdata[i][0]+'\n'))
    outputfile=open("../" + fold + "/" + term + "sum" + suf,'w')
    outputfile.writelines(l0)
    outputfile.close()
    
  
    

def main():
    '''
    T(deg):0         Z(cm):1        Xun(mm-mrad):2   Yun(mm-mrad):3   Zun(mm-mrad):4  Xn(mm-mrad):5    Yn(mm-mrad):6    
    Zn(mm-mrad):7     Xrms(mm):8       Yrms(mm):9      Zrmz(mm):10       <kE>(MeV):11           Del-Erms:12       <X>(mm):13      
    <Xpn>(mrad):14   <Y>(mm):15        <Ypn>(mrad):16  <Z>(cm):17       <Zpn>(rad):18   EZref(MV/m):19
    '''
    interestcol=[3,5,6,8]
    
    linepick(interestcol)
    
if __name__ == '__main__':
    main()
