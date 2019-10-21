import sys
import os
import numpy as np
import math

directory = r'G:\Hopkins_Laboratory\Pesticide_AcE_Inhibitor\Structures\Protonated\DFT_and_MobCal\8\Outputs'

#constants
R = 8.3144598E-3 # kJ mol^-1
T = 298.15 #kelvin

#csv containing filenames and HF energies (unscaled)
opf = open(directory+'\\Gibbs.csv', 'r')
data = opf.readlines()
opf.close()

#do not edit past this point
#curent name
cn = ''
mn = [] #master name array
me = [] #master energy array
tn = [] #temp name array
te = [] #temp energy array

for line in data[1:]:
    line = line.split(',')
    name = line[0]
    energy = float(line[1].strip())
    prefix = name[:name.rfind('_')]
    if prefix != cn: #group isomers by chemical name found in input filename
        mn.append(tn)
        me.append(te)
        tn = []
        te = []
        cn = prefix
    tn.append(name)
    te.append(energy)

#remove blank initial value in mn and me arrays
mn = mn[1:]
me = me[1:]

me2 = []
weight = []

for nrg in me:
    #get min energy
    min_e = min(nrg)
    rel_e = (np.array(nrg)-min_e)*2625.5
    #boltzmann weighting
    b = np.exp((-rel_e)/(R*T))
    pop = sum(b)
    w = b/pop
    me2.append(rel_e.tolist())
    weight.append(w.tolist())

#write filename, energy, relative energy, and botlzmann weight to file
opf = open(directory+'\\Rel_E.csv', 'w')
opf.write('Name,Energy,Relative E,Boltzmann Weight\n')
for i in range(len(mn)):
    for a in range(len(mn[i])):
        opf.write(mn[i][a]+','+str(me[i][a])+','+str(me2[i][a])+','+str(weight[i][a])+'\n')

opf.close()
