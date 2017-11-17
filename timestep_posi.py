'''
Created on Jan 29, 2015

@author: wange
This file is taking a exactly position particle information from a sort of timestepemittance.tbl file
to use this code. Run phase scan code first and got the timestep from parmela followed by various phase.
This code will generate a jitter.dat file 
Copy phase vs phase jitter of three jitter.dat files into a spreadsheet file which including str, phase, amp. Jitter study could base on these data.
In jit_phase mode TIMESTEPEMITTANCE???.TBL--->TIMESTEPEMITTANCE??????.TBL
'''
import os 
import fnmatch
path="C:\\cygwin64\\home\\wange\\UED\\1.5eleb\\par\\final\\sta\\fullrange\\str"
os.chdir(path)
if os.path.isfile('jitter.dat'):
    os.remove('jitter.dat')
target=100.0000

print os.listdir(path)
for root, dirs, files in os.walk(path):
        timestep=[elem for elem in files if fnmatch.fnmatch(elem,'TIMESTEPEMITTANCE???.TBL')]
print timestep
step=int(filter(str.isdigit,timestep[1]))-int(filter(str.isdigit,timestep[0]))
lengthts=len(timestep)
starts=int(filter(str.isdigit,timestep[0]))
print step, lengthts,starts
          

for name in timestep:
    testfile=open(name,'r')
    lines=testfile.readlines()
    testfile.close()
    steps=len(lines)-85
    print name, '  ',len(lines)
    jitterfile=open('jitter.dat','a')
    for i in range(86,len(lines)+1):
        lines[i]=lines[i].split()
        #print i,' ',  lines[i]
        if float(lines[i][1])-target==0:
            
            print str(filter(str.isdigit,name)),lines[i]
            newline=map(float,lines[i])
            newline=map("{:.10f}".format,newline)  #newline=map(lambda x:"{:.10f}".format(float(x)),newline)    or      newline=["{:.10f}".format(float(i)) for i in lines]
            # newline=map(str,newline)
            
         
            jitterfile.write(str(filter(str.isdigit,name))+'    '+'            '.join(newline)+'\n')
            #jitterfile.close()
            break
        elif float(lines[i-1][1])-target<0 and float(lines[i][1])-target>0:
            #jitterfile=open('jitter.dat','w')
            liness=map(float,lines[i-1])
            linesl=map(float,lines[i])
            #print liness, '\n',linesl
            newline=[]
            for n in range(20):
                
                newline.append(str("{:.10f}".format(((target-liness[1])/(linesl[1]-liness[1]))*(linesl[n]-liness[n])+liness[n])))
                n+=1
            jitterfile.write(str(filter(str.isdigit,name))+'    '+'            '.join(newline)+'\n')
            break
        i+=1
    jitterfile.close()
            
        
        
        
            
    


 