import sys
import os
import re

"""
enter directory containing .mout files here, wrapped with r'/file/path/here'
example:
directory = r'C:\JamFlex\MobJed\Imp_AccuracyCheck\imp100'
"""

directory = r'G:\Hopkins_Laboratory\Lipidz_Tingz\Ag_MobCalMPI_Cal\CCS_Calibration\Initial_testing'

#Please do not edit past this point
files = [x for x in os.listdir(directory) if x.lower().endswith('.mout')]

opf = open(directory+'//CCS_DATA.csv','w')
opf.write('Filename, mass, Cross Section,Standard Deviation,Number of atoms,time\n')
opf.close()

for file in files:
    opf = open(directory+'//'+file,'r')
    data = opf.read()
    opf.close()

    filename = re.findall('input file name = (.*?)\n',data)
    mass = re.findall('mass of ion = (.*?)\n',data)
    cross_section = re.findall('average TM cross section = (.*?)\n',data)
    std_dev = re.findall(r'standard deviation \(percent\) = (.*?)\n',data)
    time = re.findall(r'Job Completed in(.*?)s\n',data)
    n_atoms = re.findall(r'number of atoms =(.*?)\n',data)
    clean_up = [filename,mass,cross_section,std_dev,n_atoms,time]
    
    for parameter in range(len(clean_up)):
        if clean_up[parameter]:
            clean_up[parameter] = clean_up[parameter][-1].strip()
        else: #If these values cannot be found
            clean_up[parameter] = 'null'
    opf = open(directory+'//CCS_DATA.csv','a')
    opf.write('%s,%s,%s,%s,%s,%s\n'%(clean_up[0],clean_up[1],clean_up[2],clean_up[3],clean_up[4],clean_up[5]))
    opf.close()


    
    
