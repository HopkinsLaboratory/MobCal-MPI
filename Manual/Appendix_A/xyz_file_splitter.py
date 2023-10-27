'''
Created by CI 2023-10-11

Extractes conformers from a multi-conformers.xyz file and writes each one to its own .gjf file (Gaussian) 
'''

directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Appendix_A\CREST_Outputs'
filename = 'crest_conformers.xyz'
basename = 'Fluoxetine'

#Do not edit past this point unless you know what you are doing!

# Define file paths
import os

opf_path = os.path.join(directory, filename)
output_folder = 'Conformers'
output_path = os.path.join(directory, output_folder)
energies_path = os.path.join(output_path, 'Energies.csv')

# Read the data from the input file
with open(opf_path, 'r') as opf:
    data = opf.read()

# Creating the output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# Writing header to the energies file
with open(energies_path, 'w') as opf2:
    opf2.write('Filename,Energies,Relative Energy (kJ/mol)\n')

# Split the data into conformers
split_by = str(data.split('\n')[0].strip()) + '\n        -'
data = data.split(split_by)[1:]

# Get the minimum energy
E_min = -1 * float(data[0].split('\n')[0])

i = 1  # Set index to one for naming files

for entry in data:
    entry = entry.split('\n')
    filename_prefix = basename if not basename.endswith('_') else basename[:-1] #ensure basename ends with an underscore forreadability of resulting filenames

    # Write the .gjf file
    gjf_file_path = os.path.join(output_path, f'{filename_prefix}_{i}.gjf')
    with open(gjf_file_path, 'w') as opf3:
        opf3.write('#opt\n\n')
        opf3.write(f'{filename[:-4]}\n\n')
        opf3.write('1 1\n')  # Default charge and multiplicity. These values don't matter because you'll update them later. 
        opf3.write('\n'.join(entry[1:]))
        opf3.write('\n\n')

    # Write to the energies file
    with open(energies_path, 'a') as opf4:
        opf4.write(f'{filename_prefix}_{i}.gjf,-{entry[0]},{((-1 * float(entry[0]) - E_min) * 2625.5)}\n')

    i += 1 #update index for file naming
