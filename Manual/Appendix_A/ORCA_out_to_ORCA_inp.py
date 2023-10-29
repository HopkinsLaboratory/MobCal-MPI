# Created by CI 2022-02-19
# Converts ORCA output files (.out) to ORCA input files (.inp and .xyz for geometry)

# Where are your xyz files located?
directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Appendix_A\DFT'

# Important things for the input block
mpp = 3400          # Memory per processor in MB. Must be an integer
ncores = 8          # Number of cores to use in the calculation
charge = 1          # What is the charge?
multiplicity = 1    # What is the multiplicity?

'''Typical workflow'''
# calc_line = '! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3 CHELPG'    # DFT calcs, thermochemistry, and partial charges
calc_line = '! DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C VeryTightSCF'                            # DLPNO-CCSD(T) calcs

# optional flags
'''ESP Charges'''
ESP_Charges = False  # Do you want to calculate partial charges via ESP mapping using the ChelpG scheme?

'''Polarization'''
Polarization = False  # Do you want to calculate dipole and quadrupole moments?

'''Optional Outputs'''
write_gjf = True    # Would you like .gjf files to be made in addition to the ORCA inputs? This option is useful for visualizing files in Gaussview. 
calc_hess = False # Do you want to calculate the Hessian for the first opt step?


###########################################################################
# No touchy past this point

import os, re

# Generate a directory to write new files to
new_dir = os.path.join(directory, 'New_Inputs')
os.makedirs(new_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Generate a list of ORCA output files from the directory, excluding xyz files that are trajectories from an optimization
filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x]

# Create ORCA files
for filename in filenames:
    orca_filename = os.path.join(new_dir, f'{filename[:-4]}.inp')
    with open(orca_filename, 'w') as opf:
        # Input block
        if 'CHELPG' in calc_line and not ESP_Charges:
            ESP_Charges = True
        if 'CHELPG' not in calc_line and ESP_Charges:
            calc_line = f'{calc_line} CHELPG'

        if not calc_line.startswith('!'):
            calc_line = f'! {calc_line}'
        opf.write(f'{calc_line}\n')

        # Memory
        if not isinstance(mpp, int):
            mpp = round(mpp)
            print('The number of memory is not an integer. You cannot have fractional CPUs! Fixing it now')
        opf.write(f'%maxcore {mpp}\n\n')

        # Number of cores
        if not isinstance(ncores, int):
            ncores = round(ncores)
            print('The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
        opf.write(f'%pal nprocs {ncores} \nend\n\n')

        # ESP Charges
        if ESP_Charges:
            
            grid = 0.1 # Default is 0.3
            rmax = 3.0 # Default is 2.8

            opf.write(f'%chelpg'+'\n')
            opf.write(f'grid {grid}\n')
            opf.write(f'rmax {rmax}\n')
            opf.write('end\n\n')

        # Open shell calcs
        if multiplicity != 1:
            opf.write(f'%scf HFTyp UHF\n')
            opf.write('end\n\n')

        # Calc Hessian during optimization
        if calc_hess:
            opf.write(f'%geom\n')
            opf.write('Calc_Hess true\n')
            opf.write('end\n\n')

        # Freq cutoff
        if Polarization:
            opf.write(f'%elprop\n')
            opf.write('Dipole true\n')
            opf.write('Quadrupole true\n')
            opf.write('Polar 1\n')
            opf.write('end\n\n')

    with open(os.path.join(directory, filename), 'r') as opf:
        data = opf.read()

    # Get the geometry and write to a list called geometry
    geometry = [] # Get geometry of atoms from file
    GEOM = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data)

    # Loop through each line in the geom block, split each line, then append to the geometry list
    for i in GEOM[-1].split('\n')[2:-3]:
        i = i.split()
        if len(i) > 0:
            geometry.append(i)

    # Write the geometry to a formatted string
    fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])

    # At long last, we can write all the info to the .xyz file
    with open(orca_filename, 'a') as opf:
        opf.write(f'*xyz {charge} {multiplicity}\n')
        opf.write(fs)
        opf.write('\n*\n\n') #terminate geom block with an asterisk, followed by two blank lines

    if write_gjf:
        with open(os.path.join(new_dir, f'{filename[:-4]}.gjf'), 'w') as opf:
            opf.write('#opt\n\n')
            opf.write(f'{filename[:-4]}\n\n') #filename exclusing its extenstion 
            opf.write(f'{charge} {multiplicity}\n')
            opf.write(fs)
            opf.write('\n\n')

print('Done')
