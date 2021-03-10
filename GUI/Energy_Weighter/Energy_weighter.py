import sys
import os
import numpy as np
import math

def ew_run(self):
    csv_file = self.ui.t2le_1.text().strip()
    directory = os.path.dirname(csv_file)
    if os.path.isfile(csv_file) == False:
        self.ui.t2le_1.setStyleSheet("background-color: rgb(255, 0, 0)")
        return

    #constants
    R = 8.3144598E-3 # kJ mol^-1
    T = self.ui.t2dsb1.value()

    #csv containing filenames and HF energies (unscaled)
    opf = open(csv_file, 'r')
    data = opf.readlines()
    opf.close()

    #do not edit past this point

    #Get unique compound names
    ps = [] #prefixes
    for line in data[1:]:
        line = line.split(',')
        name = line[0]
        prefix = name[:name.rfind('_')]
        if prefix not in ps: ps.append(prefix)

    mn = [[] for x in range(len(ps))]
    me = [[] for x in range(len(ps))]

    for line in data[1:]:
        line = line.split(',')
        energy = float(line[1].strip())
        name = line[0]
        prefix = name[:name.rfind('_')]
        #Index the prefix from unique prefixes
        pi = ps.index(prefix) 
        mn[pi].append(name)
        me[pi].append(energy)

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
    opf = open(directory+'/Rel_E.csv', 'w')
    opf.write('Name,Energy,Relative E,Boltzmann Weight\n')
    for i in range(len(mn)):
        for a in range(len(mn[i])):
            opf.write(mn[i][a]+','+str(me[i][a])+','+str(me2[i][a])+','+str(weight[i][a])+'\n')

    opf.close()
