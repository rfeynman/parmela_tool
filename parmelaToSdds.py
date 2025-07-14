#!/usr/bin/env python
######################################################################
# Convert Parmela output to elegant SDDS beam format for EIC electron
#   preinjector/SRF linac start to end tracking. Note that after
#   conversion, sddsanalyzebeam can be used on the output SDDS file
#   to evaluate many beam distribution parameters.
#
# T. Satogata, July 5 2024
######################################################################
import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Constants, filenames etc
inputFile = "2024-07-05-TAPE3.TXT"
outputFile = "2024-07-05-TAPE3.sdds"
analysisFile = "2024-07-05-TAPE3-analyze.sdds"
xlsFile = "2024-07-05-TAPE3.xls"
electronMass = 510998.95 # eV/c^2
speedOfLight = 299792458 # m/s

# Read and parse input file into pandas data frame
delimiter = r'\s+'
columnNames = ['x','bgx','y','bgy','z','bgz','species','charge']
columnUnits = ['cm','--','cm','--','cm','--','-','e']
skipRows = 7
df = pd.read_csv(inputFile, sep=delimiter, header=None, names=columnNames, skiprows=skipRows)
numParticles=len(df)

# Convert x,y,z lengths to meters
df.x/=100
df.y/=100
df.z/=100

# Add new columns for total scaled momentum, xp, yp, gamma, beta, xp, yp
df['bgp'] = np.sqrt(df['bgx']**2+df['bgy']**2+df['bgz']**2)
df['p'] = df['bgp']*electronMass # eV/c
df['E'] = np.sqrt(df['p']**2+electronMass**2) # eV
df['gamma'] = df['E']/electronMass
df['beta'] = np.sqrt(1-1/df['gamma']**2)
df['xp'] = df['bgx']/df['bgp']
df['yp'] = df['bgy']/df['bgp']
df['t'] = df['z']/(df['beta']*speedOfLight) # ~33ns for 10m injector

# Remove means
xMean = df.x.mean()
yMean = df.y.mean()
zMean = df.z.mean()
tMean = df.t.mean()
df['x']-=xMean
df['y']-=yMean
df['z']-=zMean
df['dt']=df['t']-tMean

bgxMean = df.bgx.mean()*electronMass
bgyMean = df.bgy.mean()*electronMass
bgzMean = df.bgz.mean()*electronMass
ptot=np.sqrt(bgxMean**2+bgyMean**2+bgzMean**2)
ke=ptot-electronMass

print(f"means: x={1e3*xMean:.6f}mm y={1e3*yMean:.6f}mm z={zMean:.6f}m pAve={(1e-6*ptot):.3f}MeV/c keAve={(1e-6*ke):.3f}MeV")

# plt.figure(figsize=(8,5))
# plt.plot(df.x,df.bgx,marker='.',linestyle='None')
# plt.xlabel('x [m]')
# plt.ylabel('bgx [?]')
# plt.title(f'xmean={xMean:.5f}, bgxmean={bgxMean:.5f}')
# plt.grid(True)
# plt.show()
# 
# plt.figure(figsize=(8,5))
# plt.plot(df.y,df.bgy,marker='.',linestyle='None')
# plt.xlabel('y [m]')
# plt.ylabel('bgy [?]')
# plt.title(f'ymean={yMean:.5f}, bgymean={bgyMean:.5f}')
# plt.grid(True)
# plt.show()
 
with open(outputFile,'w') as OUT:
    OUT.write("SDDS5\n")
    OUT.write("&description text=\"Converted Parmela output from Erdong Wang\", contents=\"phase space\", &end\n")
    OUT.write("&parameter name=Step, description=\"Simulation step\", type=long, &end\n")
    OUT.write("&parameter name=pCentral, symbol=\"p$bcen$n\", units=\"m$be$nc\", description=\"Reference beta*gamma\", type=double, &end\n")
    OUT.write("&parameter name=Charge, units=C, description=\"Bunch charge before sampling\", type=double, &end\n")
    OUT.write("&parameter name=Particles, description=\"Number of particles before sampling\", type=long, &end\n")
    OUT.write("&parameter name=IDSlotsPerBunch, description=\"Number of particle ID slots reserved to a bunch\", type=long, &end\n")
    OUT.write("&parameter name=SVNVersion, description=\"SVN version number\", type=string, &end\n")
    OUT.write("&parameter name=SampledCharge, units=C, description=\"Sampled charge\", type=double, &end\n")
    OUT.write("&parameter name=SampledParticles, description=\"Sampled number of particles\", type=long, &end\n")
    OUT.write("&parameter name=Pass, type=long, &end\n")
    OUT.write("&parameter name=PassLength, units=m, type=double, &end\n")
    OUT.write("&parameter name=PassCentralTime, units=s, type=double, &end\n")
    OUT.write("&parameter name=ElapsedTime, units=s, type=double, &end\n")
    OUT.write("&parameter name=ElapsedCoreTime, units=s, type=double, &end\n")
    OUT.write("&parameter name=MemoryUsage, units=kB, type=long, &end\n")
    OUT.write("&parameter name=s, units=m, type=double, &end\n")
    OUT.write("&parameter name=Description, format_string=%s, type=string, &end\n")
    OUT.write("&parameter name=PreviousElementName, format_string=%s, type=string, &end\n")
    OUT.write("&column name=x, units=m, type=double,  &end\n")
    OUT.write("&column name=xp, type=double,  &end\n")
    OUT.write("&column name=y, units=m, type=double,  &end\n")
    OUT.write("&column name=yp, type=double,  &end\n")
    OUT.write("&column name=t, units=s, type=double,  &end\n")
    OUT.write("&column name=p, units=\"m$be$nc\", type=double,  &end\n")
    OUT.write("&column name=dt, units=s, type=double,  &end\n")
    OUT.write("&column name=particleID, type=ulong64,  &end\n")
    OUT.write("&data mode=ascii, &end\n")
    OUT.write("! page number 1\n")
    OUT.write("1\n") # step
    OUT.write(f"{ptot}\n") # pCentral [beta*gamma]
    OUT.write("7e-9\n") # 7 nC
    OUT.write(f"{numParticles}\n")
    OUT.write(f"{numParticles}\n")
    OUT.write("29151M\n")
    OUT.write("7e-9\n") # 7 nC
    OUT.write(f"{numParticles}\n")
    OUT.write("0\n") # Pass
    OUT.write(f"{zMean}\n") # PassLength [m]
    OUT.write("0.000000000000000e+00\n") # PassCentralTime [s]
    OUT.write("0.000000000000000e+00\n") # ElapsedTime [s]
    OUT.write("0.000000000000000e+00\n") # ElapsedCoreTime [s]
    OUT.write("0\n") # MemoryUsage [kB]
    OUT.write("0.000000000000000e+00\n") # s [m]
    OUT.write("\"\"\n") # Description
    OUT.write("CHARGE\n") # previousElementName
    OUT.write(f"{numParticles}\n")
    # x1 xp1 y1 yp1 t1 bgp1 dt1 particleID1
    for index, row in df.iterrows():
        OUT.write(f"{df.at[index,'x']} {df.at[index,'xp']} {df.at[index,'y']} {df.at[index,'yp']} {df.at[index,'t']} {df.at[index,'bgp']} {df.at[index,'dt']} {index+1}\n")

os.system(f"sddsanalyzebeam {outputFile} {analysisFile}")
os.system(f"sdds2spreadsheet {analysisFile} {xlsFile} -excel")
