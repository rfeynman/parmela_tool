'''
Created on Jul 19, 2017
This code for compbine the table. TBL file from parmela
@author: Erdong Wang
combine multiple colume into sigle TIMESTEP file. choose the coleum and phase steps and run the code
'''

import numpy as np
from numpy import genfromtxt, dtype, loadtxt
import string



fold = "560"
term = "TIMESTEPEMITTANCE"
suf = ".TBL"
listrange = list(range(50, 110, 10))
dataline=84


def linepick(pick):
    sumdata=[]
    datanum=[]
    for intr in listrange[1:]:
        with open("../" + fold + "/" + term + str(intr) + suf) as f:
            l = [line.split() for line in f]
        f.close()
        usefuldata = l[dataline:]

        # content=genfromtxt("../"+fold+"/"+term+"90"+suf,delimiter='/t',names=True,dtype=None)
        headname = [elem for elem in usefuldata[0]][1:]
        usefuldata[0]=headname
        col_usefuldata={usefuldata[0][i]:[line[i] for line in usefuldata] for i in range(len(usefuldata[0]))}
        for i in range(len(headname)):
            col_usefuldata[headname[i]][0]=headname[i]+str(intr)
            
        for elem in pick:
            sumdata.append(col_usefuldata[headname[elem]])
        #print(col_usefuldata[headname[0]])
        datanum.append(len(usefuldata))
    minnum=min(datanum)
    with open("../" + fold + "/" + term + str(listrange[0]) + suf) as f0:
        l0=f0.readlines()
    #print(l0[84])
    f0.close()
    
    #print(l0[84].strip('\n'),str(sumdata[0][0]),"\n",'    '.join([l0[84].strip('\n'),str(sumdata[0][0])]))

    
    for k in range(0,minnum):
        for i in range(len(sumdata)):
            
            l0[k]='    '.join([l0[k+dataline].strip('\r\n'),str(sumdata[i][k])+'\r\n'])
   
    del l0[minnum+dataline:]
   
    inslabel=l0.index('Ezref(MV/m)\r\n')
	

    for i in range(len(sumdata)):
        l0.insert(inslabel+i+1,str(sumdata[i][0]+'\r\n'))
        l0.insert(inslabel+2*i+25,str(sumdata[i][0]+'\r\n'))
        l0.insert(inslabel+3*i+49,str(sumdata[i][0]+'\r\n'))
    outputfile=open("../" + fold + "/" + term + "sum" + suf,'w')
    outputfile.writelines(l0)
    outputfile.close()
    
  
    

def main():
    interestcol=[10,19]
    
    linepick(interestcol)
    
if __name__ == '__main__':
    main()
