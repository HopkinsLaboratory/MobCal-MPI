'''
Program to abstract the calculation of spin multiplicity for target chemical.
'''

import sys
import os
import ast

# Update directory from CLI
if (len(sys.argv) > 1):
    #Get passed in charge and file
    ionFile = sys.argv[1]
    ionCharge = ionFile.split('.')[0]
else:
    print("A file input is required...", file=sys.stderr)
    sys.exit(1)

# Ion Multiplicity calculations are difficult and need further work...
ionMultiplicity = 1
ionMultiplicityFile="spin_multiplicity.config"

#Create blank dictionary
ionDictionary = {}

#If file exists
if(os.path.exists(ionMultiplicityFile)):
    #Read data from file and update ion multiplicity value
    with open(ionMultiplicityFile, "r") as data:
        ionDictionary = ast.literal_eval(data.read())
        #If the ion charge key is in the dictionary
        if ionCharge in ionDictionary.keys():
            #Update ion multiplicity value
            ionMultiplicity = ionDictionary[ionCharge]

#Place ion multiplicity value into dictionary
ionDictionary[ionCharge] = ionMultiplicity

#Write updated dictionary to file
with open(ionMultiplicityFile, "w") as data:
    data.write(str(ionDictionary))

print("If run fails, update value '" + ionCharge + "' in '" + ionMultiplicityFile + "' for future runs.", file=sys.stderr)

#Print result
print(ionMultiplicity)
