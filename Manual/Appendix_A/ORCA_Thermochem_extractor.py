# Created by CI 2021-04-15
# Extracts thermochemical quantities directly from ORCA output files

import os
import re
import sys

# User-defined settings
directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Appendix_A\DFT'

# No touchy past this point unless you know what you are doing. 

thermo_properties = [
    'Filename',
    'Imaginary Frequencies',
    'Electronic energy',
    'ZPE correction',
    'Thermal correction',
    'Enthalpy correction',
    'Gibbs Correction',
    'Total ZPE',
    'Total Thermal Energy',
    'Total Enthalpy',
    'Total Entropy',
    'Total Gibbs Energy'
]

# Format the header for consistent spacing 
header = '{}\n'.format(','.join(['{:<25}'] * len(thermo_properties)))

#Create output file and write header to it
output_csv = os.path.join(directory, 'Thermo_data.csv')

try:
    with open(output_csv, 'w') as opf:
        opf.write(header.format(*thermo_properties))

except PermissionError:
    print('A file with the same name is already open. Close it and rerun the code')
    sys.exit()

# Function to extract a numerical property from the ORCA output file
def extract_property(data, pattern):
    try:
        result = float(re.findall(pattern, data)[-1].strip())
    except:
        result = 12345.0  # Placeholder for missing values
    return result

for filename in [x for x in os.listdir(directory) if x.lower().endswith('.out')]:
    with open(os.path.join(directory, filename), 'r') as opf:
        data = opf.read()

    # Extract thermochemical properties from the ORCA output file
    E_el = extract_property(data, 'FINAL SINGLE POINT ENERGY(.*?)\n')
    ZPE_corr = extract_property(data, 'Zero point energy                ...(.*?)Eh')
    Thermal_corr = extract_property(data, 'Total thermal correction(.*?)Eh')
    H_corr = extract_property(data, 'Thermal Enthalpy correction       ...(.*?)Eh')
    S_tot = extract_property(data, 'Final entropy term                ...(.*?)Eh')
    G_corr = extract_property(data, r'G-E\(el\)                           ...(.*?)Eh')
    E_thermal = extract_property(data, 'Total thermal energy(.*?)Eh')
    E_Enthalpy = extract_property(data, 'Total Enthalpy                    ...(.*?)Eh')
    E_Gibbs = extract_property(data, 'Final Gibbs free energy         ...(.*?)Eh')
    E_ZPE = E_el + ZPE_corr
    n_imag = data.count('***imaginary mode***')

    if 12345.0 in [E_el, ZPE_corr, Thermal_corr, H_corr, G_corr, E_ZPE, E_thermal, E_Enthalpy, S_tot, E_Gibbs]:
        print(f'{filename} is missing thermochemistry. Writing -12345 as a placeholder for missing value\n')

    if n_imag > 0:
        print(f'{filename} contains {n_imag} imaginary frequencies.\n')

    # Prepare the values to be written
    values = [filename, n_imag, E_el, ZPE_corr, Thermal_corr, H_corr, G_corr, E_ZPE, E_thermal, E_Enthalpy, S_tot, E_Gibbs]

    # Create a format string for consistent spacing
    format_str = '{}'.format('{:<25},' * (len(values) - 1) + '{:<25}\n')

    # Append the extracted properties to the CSV file
    with open(output_csv, 'a') as opf:
        opf.write(format_str.format(*values))

print('Done')
