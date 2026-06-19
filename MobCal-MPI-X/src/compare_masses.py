#Script adapted from Featherstone GUI functions
import os
import shutil
import sys
import numpy as np

from numpy import dot
from numpy.linalg import norm

from rearrange_molecules import GEOM
from connectivity import find_connectivity

#Declare user input variables to update from CLI
directory = ""
threshold = 0

# Update directory from CLI
if (len(sys.argv) > 1):
    directory = sys.argv[1]

# Update similarity threshold from CLI
if (len(sys.argv) > 2):
    threshold = float(sys.argv[2])
else:
    print("Directory and similarity threshold inputs are required...")
    sys.exit(1)

def determine_similarity(directory, similarity_threshold):
    #Get all GJF files in target directory
    Files = sorted([x for x in os.listdir(directory) if x.endswith('.gjf')])

    mwl = [] #mass weighted coordinate list
    for File in Files:
        filePath = f'{directory}/{File}'
        try:
            with open(filePath,'r') as opf:
                lines = opf.readlines()
                geometry = []
                labels = []
                for line in lines:
                    line = line.split()
                    if len(line) != 4: continue
                    try: 
                        geometry.append([float(line[1]),float(line[2]),float(line[3])])
                        labels.append(line[0])
                    except ValueError: continue
                geometry = np.array(geometry)
                mw = find_connectivity(labels,geometry)
                mwl.append(mw)
        except KeyError:
            print("File error resulted during connectivity step. Skipping: " + File, file=sys.stderr)
            # If the file exists, remove it to prevent cascading errors
            if(os.path.exists(filePath)):
                os.remove(filePath)
            # Remove from list of files
            Files.remove(File)
        except FileNotFoundError:
            try:
                with open(f'{directory}/Error.txt','a+') as erropf:
                    erropf.writelines(f'{File} \n')
            except FileNotFoundError:
                with open(f'{directory}/Error.txt','w') as erropf:
                    erropf.writelines("Unable to find the following file(s), please double check if they exist. \n")
                    erropf.writelines(f'{File} \n')
            
                
    failed = []
    for index1, mw1 in enumerate(mwl):
        if index1 in failed: continue
        for index2, mw2 in enumerate(mwl):
            if index2 <= index1 or index2 in failed: continue
            #similarity = round(dot(mw1, mw2)/(norm(mw1)*norm(mw2)),8)
            try: similarity = round(abs(1-np.average(np.array(mw1)/np.array(mw2))),8)
            except ValueError: continue #If structures of two different sizes are compared
            if similarity < (100-similarity_threshold)/100:
                #print(Files[index1],Files[index2],mw1,mw2,similarity)
                failed.append(index2)

    for index in sorted(failed,reverse=True):
        del Files[index]

    try: os.mkdir(f'{directory}/uniques')
    except FileExistsError: pass

    for File in Files:
        shutil.copy(f'{directory}/{File}',f'{directory}/uniques/{File}')

#Call function to determine overall similarity given inputs
determine_similarity(directory, threshold)

