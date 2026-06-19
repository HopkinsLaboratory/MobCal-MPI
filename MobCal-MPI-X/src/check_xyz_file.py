import sys
import math

errorNumber = 0

def ReadXYZFile(filename):
    #Create file handle
    fileHandler = open(filename, "r")

    #Read all contents and split
    content = fileHandler.readlines()

    #Close file handle
    fileHandler.close()

    for i in range(0, len(content)):
        content[i] = content[i].replace('\n', '')
    
    return content

if (len(sys.argv) > 1):
    #Open chemical and solvent files
    chemical = ReadXYZFile(sys.argv[1])
    lineOne = []
    lineTwo = []
    lineThree = []
    firstLinePass = False
    secondLinePass = True
    thirdLinePass = False
    
    #Try to access all elements of the file
    if len(chemical) > 3:
        lineOne = chemical[0].split()
        lineTwo = chemical[1].split()
        lineThree = chemical[2].split()
    else:
        errorNumber = 2

    #If there is nothing in line one or the value does not equal the total number of atoms
    if (errorNumber == 0) and ((len(lineOne) != 1) or (not lineOne[0].isdigit()) or (float(lineOne[0]) != (len(chemical) - 2))):
        errorNumber = 3

    #If there is nothing in line one or the value does not equal the total number of atoms
    if (errorNumber == 0) and ((len(lineThree) == 0) or (not lineThree[0].isalpha())):
        errorNumber = 3
else:
    errorNumber = 1

#If all checks pass
if errorNumber == 0:
    #Print out XYZ to be processed by other scripts
    print("XYZ")
#Otherwise, print to the error stream
elif errorNumber == 1:
    print("ERROR: Input file not provided to checker.", file=sys.stderr)
elif errorNumber == 2:
    print("ERROR: Input file did not have enough lines to be format.", file=sys.stderr)
elif errorNumber == 3:
    print("ERROR: Input file format did not match.", file=sys.stderr)
