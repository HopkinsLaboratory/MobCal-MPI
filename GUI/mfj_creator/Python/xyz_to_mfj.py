import os
import random
from PyQt6.QtWidgets import QMessageBox

def xyz_to_mfj(self, directory, xyz, key, mfj, charge, parameters):
	
	def write_error(filename, atom_number):
		error_file_path = os.path.join(mfj[:mfj.rfind('\\')+1], 'Errors.csv')
		mode = 'w' if not os.path.isfile(error_file_path) else 'a'
		
		with open(error_file_path, mode) as error:
			error.write(f'{filename},{atom_number}\n')

	# Read xyz file
	with open(xyz, 'r') as xyz_file:
		xyz_lines = xyz_file.readlines()

	# Read key file
	with open(key, 'r') as key_file:
		key_lines = key_file.readlines()

	# Read mass.prm file
	with open(os.path.join(directory, 'mass.prm'), 'r') as mass_file:
		mass_lines = mass_file.readlines()

	# Read vdw.prm file
	with open(os.path.join(directory, 'vdw.prm'), 'r') as vdw_file:
		vdw_lines = vdw_file.readlines()

	atom_num = int(xyz_lines[0].split()[0])
	xyz_data = [line.split() for line in xyz_lines[1:atom_num + 1]]
	key_data = [line.split() for line in key_lines[:atom_num]]

	#dictionaries for storing information from .prm files
	atom_info = {}
	mass_info = {}
	for line in mass_lines:
		parsed = line.split()
		atom_info[parsed[1]] = parsed[2]
		mass_info[parsed[2]] = parsed[6]

	vdw_info = {}
	for line in vdw_lines:
		parsed = line.split()
		vdw_info[parsed[0]] = [parsed[1], parsed[2], parsed[3], parsed[4]]

	seed = random.randint(1000000, 1000000000)
	seed = -seed #seed needs to be negative for RANLUX to initialize

	#Create .mfj file
	with open(mfj, 'w') as output_file:
		#write file header
		output_file.write('{}\n'.format(mfj.split('.')[0].replace('/', '\\').split('\\')[-1]))
		output_file.write(
			'1\n'
			f'{atom_num}\n'
			'ang\n'
			f'{charge}\n'
			f'{parameters[4]}\n'
			f'{parameters[0]} '
			f'{parameters[1]} '
			f'{parameters[2]} '
			f'{parameters[3]} '
			f'{seed} '
			f'{parameters[5]}\n'
		)

		#get atomic mass, and vdW paramaters from parameter files
		for i in range(atom_num):
			atom_mass = mass_info.get(atom_info.get(xyz_data[i][5], ''),'')
			vwd_w = vdw_info.get(atom_info.get(xyz_data[i][5], ''), ['', '', '', ''])

			#if parameter is missing, write relevant information to Errors.csv file
			if not atom_mass:
				try:
					write_error(mfj.split('.')[0].replace('/', '\\\\').split('\\\\')[-1], i + 1)
					print('No mass and vdw for atom label: {} in file: {}\n'.format(i + 1, mfj.split('.')[0].replace('/', '\\').split('\\')[-1]))

				except PermissionError:
					self.error_popup('critical','sdf2tinkerxyz Error',f'Cannot write errors in assigning in MM2 atom types to the Errors.csv file, as it is open in another window. Please close the Errors.csv file and rerun the mfj creation module.')
					return			

			# write atomic number, xyz coordinates, atomic mass, and vdW parameters to .mfj file 
			output_file.write(
				f'{xyz_data[i][2]:>10}   {xyz_data[i][3]:>10}   {xyz_data[i][4]:>10}   '
				f'{atom_mass:>7}   {key_data[i][2]:>10}   '
				f'{vwd_w[0]:>5}   {vwd_w[1]:>5}   {vwd_w[2]:>5}   {vwd_w[3]:>5}\n'
			)



