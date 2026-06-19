# Modified from MobCal-MPI repository
# Converts Gaussian input (.gjf) files to ORCA input files

import sys

# Where are your xyz files located?
directory = ""
mpp = 3200  # Memory per processor in MB. Must be an integer
ncores = 8  # Number of cores to use in the calculation
charge = 1  # What is the charge?
multiplicity = 1  # What is the multiplicity?

# Update directory from CLI
if (len(sys.argv) > 1):
    directory = sys.argv[1]

# Update memory per processor from CLI
if (len(sys.argv) > 2):
    mpp = int(sys.argv[2])

# Update number of cores from CLI
if (len(sys.argv) > 3):
    ncores = int(sys.argv[3])

# Update charge from CLI
if (len(sys.argv) > 4):
    charge = int(sys.argv[4])
    
# Update multiplicity from CLI
if (len(sys.argv) > 5):
    multiplicity = int(sys.argv[5])
else:
    print("Directory, memory per processor, number of cores, charge, and multiplicity are required inputs...")
    sys.exit(1)

# Important things for the input block
calc_line = '! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3 CHELPG'  # Must be a string
# calc_line = '! DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightSCF  '  # Must be a string

# ESP Charges
ESP_Charges = True  # Do you want custom parameters in the ChelpG scheme?

# Ask orca to call an external xyz file to run the calculation (True), or write the xyz coordinaties to the .inp file (False)
write_xyz = False

# No touchy past this point!!!

# ESP grid parameters
grid = 0.1  # Default is 0.3
rmax = 3.0  # Default is 2.8

import os

# Generate a directory to write new files to
new_dir = os.path.join(directory, 'ORCAInputs')
os.makedirs(new_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Generate a list of .gjf files from the directory
filenames = [x for x in os.listdir(directory) if x.lower().endswith('.gjf')]

# Create ORCA files
for filename in filenames:
    orca_filename = os.path.join(new_dir, f'{filename[:-4]}.inp')
    with open(orca_filename, 'w') as opf:
        # Input block
        if not calc_line.startswith('!'):
            print(f'The input line in {filename} does not start with a ! Adding it now')
            calc_line = '! ' + calc_line
        opf.write(f'{calc_line}\n')

        # Memory
        if not isinstance(mpp, int):
            print(f'The number of memory is not an integer. You cannot have fractional CPUs! Fixing it now')
            mpp = round(mpp)
        opf.write(f'%maxcore {mpp}\n\n')

        # Number of cores
        if not isinstance(ncores, int):
            print(f'The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
            ncores = round(ncores)
        opf.write(f'%pal nprocs {ncores} \nend\n\n')

        # ESP Charges
        if ESP_Charges:
            opf.write(f'%chelpg\n')
            opf.write(f'grid {grid}\n')
            opf.write(f'rmax {rmax}\n')
            opf.write('end\n\n')

    with open(os.path.join(directory, filename), 'r') as opf:
        lines = opf.readlines()

    geometry = []  # Get geometry of atoms from file
    for line_index, line in enumerate(lines):
        split_line = line.split()
        if len(split_line) == 4 and all('.' in x for x in split_line[1:4]):
            if not geometry:
                geom_index = line_index
            geometry.append(line.split())

    fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])

    if write_xyz:
        with open(os.path.join(new_dir, f'{filename[:-4]}.xyz'), 'w') as opf:
            opf.write(f'{len(geometry)}\n\n{fs}\n\n')
    else:
        with open(orca_filename, 'a') as opf:
            opf.write(f'*xyz {charge} {multiplicity}\n{fs}\n*\n\n')

