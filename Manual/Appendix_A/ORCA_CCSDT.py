# Created by CI 2021-04-20
#Extracts CCSD(T) energies from ORCA output files

import os
import re
import sys

# Define the directory where ORCA output files are located
directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Appendix_A\DLPNO_CCSDT'

# No touchy past this point

properties = [
    'Filename',
    'DLPNO-CCSD(T) energy',
]

# Format the header for consistent spacing 
header = '{}\n'.format(','.join(['{:<25}'] * len(properties)))

#Create output file and write header to it
output_csv = os.path.join(directory, 'DLPNO_CCSDT_energies.csv')

try:
    with open(output_csv, 'w') as opf:
        opf.write(header.format(*properties))

except PermissionError:
    print('A file with the same name is already open. Close it and rerun the code')
    sys.exit()

# Function to extract a numerical property from the ORCA output file
def extract_CCSDT(data, pattern):
    try:
        E_CCSDT = float(re.findall(pattern, data)[-1].strip())
    
    except Exception as e:
        print(f'{filename} has an error: {e}. Writing 12345.0 as a placeholder for the missing value.')
        E_CCSDT = 12345.0
        
    return E_CCSDT

for filename in [x for x in os.listdir(directory) if x.lower().endswith('.out')]:
    with open(os.path.join(directory, filename), 'r') as opf:
        data = opf.read()

    # Extract thermochemical properties from the ORCA output file
    CCSDT = extract_CCSDT(data, r'E\(CCSD\(T\)\)                                 ...(.*?)\n')

    # Prepare the values to be written
    values = [filename, CCSDT]

    # Create a format string for consistent spacing
    format_str = '{}'.format('{:<25},' * (len(values) - 1) + '{:<25}\n')

    # Append the extracted properties to the CSV file
    with open(output_csv, 'a') as opf:
        opf.write(format_str.format(*values))

print('Done')
