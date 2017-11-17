#!/usr/bin/env python
import os, sys, string, pdb
start_input=100
end_input=200
intv=5
for i in range(start_input, end_input+1,intv):
	os.system('echo '+str(i)+'>>temp')
	os.system('sed -n $,+0p TIMESTEPEMITTANCE'+str(i)+'.TBL >>temp')