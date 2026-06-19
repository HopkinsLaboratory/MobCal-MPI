import sys
import math

def ReadXYZFile(filename):
    #Create file handle
    fileHandler = open(filename, "r")

    #Read all contents and split
    content = fileHandler.read().split()

    #Close file handle
    fileHandler.close()

    #Convert all numerical data and determine atom spacing
    characterIndices = []
    for i in range(0, len(content)):
        try:
            content[i] = float(content[i])
        except:
            content[i] = content[i]
            characterIndices = characterIndices + [i]

    characterIndices = characterIndices + [len(content)]

    data = [[content[0]]]

    for i in range(0, len(characterIndices)-1):
        data = data + [content[characterIndices[i]:characterIndices[i+1]]]

    return data

def WriteXYZFile(filename, data):

    handle = open(filename, "w")

    for i in range(0, len(data)):            
        for j in range(0, len(data[i])):
            handle.write(str(data[i][j]) + "    ")
            
        handle.write("\n")
        if(len(data[i]) == 1):
            handle.write("\n")
        elif( i < (len(data)-1) and len(data[i+1]) == 1):
            handle.write("\n")
           
    handle.close()
    
    return

def CalculateRadius(atom):
    return math.sqrt(pow(atom[1],2) + pow(atom[2],2) + pow(atom[3],2))

if (len(sys.argv) > 3):
    #Get all inputs
    outputFile = sys.argv[3]

    #Open chemical and solvent files
    chemical = ReadXYZFile(sys.argv[1])
    solvent = ReadXYZFile(sys.argv[2])

    #Hold calculated radius
    tempRadius = 0
    
    #Determine the maximum radius value
    maxChemRadius = 0
    for i in range(1, len(chemical)):
        tempRadius = CalculateRadius(chemical[i])
        if (tempRadius > maxChemRadius):
            maxChemRadius = tempRadius

    maxSolRadius = 0
    for i in range(1, len(solvent)):
        tempRadius = CalculateRadius(solvent[i])
        if (tempRadius > maxSolRadius):
            maxSolRadius = tempRadius

    #Double radius to account for potential overlap
    maxChemRadius += 1 + maxSolRadius
    
    #Shift all solvent 
    for i in range(1, len(solvent)):
        if(len(solvent[i]) > 1):
            solvent[i][1] += maxChemRadius

    #Update ion count
    chemical[0][0] = int(chemical[0][0] + solvent[0][0])
    dimer = chemical + solvent[1:len(solvent)]

    WriteXYZFile(outputFile, dimer)
