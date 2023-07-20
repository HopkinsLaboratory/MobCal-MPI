import os
import random

def xyz_to_mfj(directory,xyz,key,mfj,charge,parameters):
	#Get the filename
	filename = (mfj.split("."))[0].replace('/','\\')
	filename = filename[filename.rfind('\\')+1:]
	
	#open files and grab all data
	opf = open(xyz,"r")
	xyz_lines = opf.readlines()
	opf.close()

	opf = open(key,"r")
	key_lines = opf.readlines()
	opf.close()

	opf = open(directory+"mass.prm","r")
	mass_lines = opf.readlines()
	opf.close()

	opf = open(directory+"vdw.prm","r")
	vdw_lines = opf.readlines()
	opf.close()

	#Read the xyz file and obtain the number of atoms and get all x,y,z data
	atom_num = int((xyz_lines[0].split()[0]))
	xyz_data = []
	for i in range(1,1+atom_num):
		xyz_data.append(xyz_lines[i].split())		

	#Get the atom info in the key file
	key_data = []
	for i in range(0,atom_num):
		key_data.append(key_lines[i].split())

	#Get the atom and mass info from the mass file
	atom_info = {}
	mass_info = {}
	for i in range(len(mass_lines)):
		parsed = mass_lines[i].split()
		atom_info[(parsed)[1]] = parsed[2]
		mass_info[(parsed)[2]] = parsed[6]

	#Get all relevant vdw info from vdw file
	vdw_info = {}
	for i in range(len(vdw_lines)):
		parsed = vdw_lines[i].split()
		vdw_info[(parsed)[0]] = [parsed[1],parsed[2],parsed[3],parsed[4]]

	opf = open(mfj,"w")
	seed = random.randint(1000000,1000000000)
	seed = -(seed)
	spacing = '%s\n%s\n%s\n%s\n%s\n%s\n%s %s %s %s %s %s\n' #regex for header
	#We want to write the filename, 1, number of atoms, ang, calc, and 1 to the header
	opf.write(spacing%(filename,'1',str(atom_num),'ang',charge,parameters[4],parameters[0],parameters[1],parameters[2],parameters[3],seed,parameters[5])) # changed in V2
	for i in range(atom_num):
		spacing = '%10s	   %10s	   %10s	   %7s	  %10s	  %5s	 %5s	%5s	   %5s\n' #regex for line
		try: #Get the mass for the atom
			atom_mass = str(mass_info[atom_info[xyz_data[i][5]]])
		except (NameError,KeyError): #If the atom type cannot be determined in the mass file throw error
			atom_mass = ''
			if os.path.isfile(mfj[:mfj.rfind('\\')+1]+'\Errors.csv') == False:
				error = open(mfj[:mfj.rfind('\\')+1]+'\Errors.csv','w')
				error.write('Filename,Atom Label\n')
				var = filename+','+str(i+1)+'\n'
				error.write(var)
				error.close()
			else:
				try:
					error = open(mfj[:mfj.rfind('\\')+1]+'\Errors.csv','a')
					var = filename+','+str(i+1)+'\n'
					error.write(var)
					error.close()
				except PermissionError:
					print('Cannot write to error file, please close it!')
			print('No mass and vdw for atom label: '+str(i+1)+' in file: '+filename+'\n')
		try: #Get the van der Waals forces for the atom
			vwd_w = [vdw_info[atom_info[xyz_data[i][5]]][0],vdw_info[atom_info[xyz_data[i][5]]][1],vdw_info[atom_info[xyz_data[i][5]]][2],vdw_info[atom_info[xyz_data[i][5]]][3]]
		except (NameError,KeyError): #If the atom type cannot be determined in the vdw file throw error
			vwd_w = ['','','','']
		#We want to write the relevant atom data to the file being: x, y, z, atom mass, key data, vdw force 1, vdw force 2, vdw force 3, vdw force 4
		opf.write(spacing%(xyz_data[i][2],xyz_data[i][3],xyz_data[i][4],atom_mass,key_data[i][2],vwd_w[0],vwd_w[1],vwd_w[2],vwd_w[3]))
	opf.close()
	
