# Created by CI 2021-04-15
# Extracts thermochemical quantities directly from ORCA output files

import os
import re
import sys
import numpy as np
import pandas as pd

# User-defined settings
directory = r'G:\Hopkins_Laboratory\Protonation_Induced_Chirality_v2\Verapamil_Impurity_A_z2\SR\DFT\1mer'

# No touchy past this point

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
    'Total Gibbs Energy',
    'Relative Gibbs Energy'
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

Gibbs_list = []

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

    Gibbs_list.append(E_Gibbs)

    if 12345.0 in [E_el, ZPE_corr, Thermal_corr, H_corr, G_corr, E_ZPE, E_thermal, E_Enthalpy, S_tot, E_Gibbs]:
        print(f'{filename} is missing thermochemistry. Writing 12345 as a placeholder for missing value')

    if n_imag > 0:
        print(f'{filename} contains {n_imag} imaginary frequencies.')

    # Prepare the values to be written
    values = [filename, n_imag, E_el, ZPE_corr, Thermal_corr, H_corr, G_corr, E_ZPE, E_thermal, E_Enthalpy, S_tot, E_Gibbs]

    # Create a format string for consistent spacing
    format_str = '{}'.format('{:<25},' * (len(values) - 1) + '{:<25}\n')

    # Append the extracted properties to the CSV file
    with open(output_csv, 'a') as opf:
        opf.write(format_str.format(*values))


# Calculate the minimum Gibbs energy
min_Gibbs = np.min(Gibbs_list)

# Read the CSV into a pandas DataFrame
df = pd.read_csv(output_csv)

# Check for the correct column name in the DataFrame
gibbs_column = next((col for col in df.columns if 'Total Gibbs Energy' in col), None)

if gibbs_column is not None:
    # Calculate relative energy column and add it to the DataFrame
    df['Relative Gibbs Energy'] = (df[gibbs_column] - np.min(Gibbs_list)) * 2625.5

    # Write the updated DataFrame back to the CSV
    df.to_csv(output_csv, index=False)
    print('Done')
else:
    print('Column "Total Gibbs Energy" not found in DataFrame.')
    df['Relative Gibbs Energy'] = 12345.0



