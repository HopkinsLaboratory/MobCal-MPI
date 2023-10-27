import os
import re
import subprocess
import time
from shutil import copyfile
from mfj_creator.Python.xyz_to_mfj import *
from PyQt5.QtWidgets import QMessageBox

def run(directory, csv, sdf2xyz2sdf_Directory, charge, parameters):
	start_time = time.time()

	def get_files(file):
		return file.lower().endswith(('.log', '.out')) and '_atom' not in file.lower() #_atom needed for ORCA jobs that use pseudopotentials add _atom(atomic number) to the filename, which are not useful for .mfj conversion

	def estimated_runtime(number_of_files):
		#about 1.5 seconds per file, printed out in converted to minutes
		return round((number_of_files * 1.5) / (60), 2)
	
	def error_popup(urgency, type, message):
		'''Sends a popup tot the GUI if proc'd. Usage is urgency level, error title, and message text'''
		msg = QMessageBox()
		msg.setWindowTitle(type)
		msg.setText(message)

		if urgency == 'warning' or 'Warning':
			msg.setIcon(QMessageBox.Warning)
		else:
			msg.setIcon(QMessageBox.Critical)

		msg.exec_()

	logs = []

	if csv != '':
		with open(csv, 'r') as opf:
			logs = [x for x in opf.read().split('\n') if get_files(x)]
	else:
		logs = [x for x in os.listdir(directory) if get_files(x)]

	print(f'Process estimated to take {estimated_runtime(len(logs))} minutes.')

	directory = directory.rstrip('\\') + '\\'

	#make MobCal_inputs directory (unless it already exists)
	try:
		os.mkdir(directory + 'Mobcal_Inputs')
	except FileExistsError:
		#if directory already exists, delete all existing contents
		[os.remove(os.path.join(os.path.join(directory, 'Mobcal_Inputs'), file)) for file in os.listdir(os.path.join(directory, 'Mobcal_Inputs'))]
	
	# Check if both Gaussian and ORCA files are present.
	gaussian_logs = any(log.lower().endswith('.log') for log in logs)
	orca_logs = any(log.lower().endswith('.out') for log in logs)

	if gaussian_logs and orca_logs:
		error_popup('critical','File Error','Gaussian files (.log) and ORCA files (.out) are both present in the specified directory. Please put the Gaussian and ORCA files into separate directories and re-run the mfj creator.')
		return
	elif not gaussian_logs and not orca_logs:
		error_popup('critical','File Error','There are no Gaussian files (.log) or ORCA files (.out) in the specified directory')
		return

	# Copy files to a new directory to avoid changing the originals
	for file in logs:
		try:
			copyfile(os.path.join(directory, file), os.path.join(directory, 'Mobcal_Inputs', file))
		except PermissionError: 
			pass
	
	directory = os.path.join(directory, 'Mobcal_Inputs')

	# Check if Open Babel 2.4.1 is installed
	try:
		version_info = subprocess.check_output(['babel', '-V'], stderr=subprocess.STDOUT, text=True)
		version_lines = version_info.splitlines()
		openbabel_version = version_lines[0].split()[2].strip() if version_lines[0].startswith('Open Babel') else None

		if openbabel_version not in ['2.4.1', '2.3.90']: #2.3.90 is an exception because the ORCA version of Avogadro preinstalls this versino of babel.  
			error_popup('critical','OpenBabel Error',f'Open Babel version: {openbabel_version} was found on your PC. Please uninstall it, and install v2.4.1 as recommended in the manual.')
			return
		
	except subprocess.CalledProcessError as e:
		print(f'Subprocess error when calling babel: {e.output}')
		return

	except FileNotFoundError:
		error_popup('critical','OpenBabel Error','Open Babel could not be found on your PC. Please install v2.4.1 as stated in the manual and ensure the babel.exe is available in your system\'s PATH')
		return

	#check for missing python files
	required_python_files = os.path.join(os.getcwd(), 'mfj_creator', 'Python')
	required_files = ['xyz_to_mfj.py', 'mass.prm', 'vdw.prm']
	missing_files = []

	for file in required_files:
		file_path = os.path.join(required_python_files, file)
		if not os.path.isfile(file_path): 
			missing_files.append(file)

	if missing_files: #checks if list is populated. If populated, will return true. If empty, will return false. 
		missing_files_text = ', '.join(missing_files)
		error_popup('critical','File Dependency Error',f'The following required Python files are missing from {required_python_files} : {missing_files_text}. Please redownload from the GitHub Repo and do not remove anything.')
		return

	print(f'If any errors are encountered, they will be written to: {directory}Errors.csv')
	if os.path.isfile(os.path.join(directory, 'Errors.csv')): 
		os.remove(os.path.join(directory, 'Errors.csv')) #remove errors.csv if it already exists

	ESP = [] #empty list for storing ESP charges from each output file

	# Workflow for handing Gaussian files: Extract ESP data, convert .log to .sdf
	if gaussian_logs:
		
		print('Extracting ESP info from Gaussian logs.\n')
		# Extracting ESP data from Gaussian .log files
		for file in logs:
			with open(os.path.join(directory, file), 'r') as opf:
				data = opf.read().replace('\n', 'ggez')

			try:
				data = re.findall(r'ESP charges:(.*?)Sum of ESP charges', data)[-1]
				data = data.split('ggez')
				data = [x.split()[2] for x in data if len(x.split()) == 3]
				ESP.append(data)
			except:
				error_popup('critical','ESP Charge Error',f'{file} is missing ESP data, did it finish correctly?\nGUI will now exit.')
				return

		print('ESP succesfully extracted from Gaussian logs.\n')

		print('Converting Gaussian .log files to .sdf via Open Babel\n')
		babel_i = os.path.join(directory, '*.log')
		babel_o = os.path.join(directory, '*.sdf')
		command = f'babel "{babel_i}" -osdf "{babel_o}"'

		try:
			convert_babel = subprocess.check_output(command, shell=True)
		except subprocess.CalledProcessError:
			try:
				print(f'Subproccess call to OpenBabel failed. Attempting to execute {command} via os.system')
				os.system(command) #alternative to subproccess call in case it fails.		
			except:
				error_popup('critical','OpenBabel Conversion Error','Encountered an error using Babel to convert the Gaussian .log files to .sdf files. Is Babel installed?')
				return

		print('Gaussian .log files were successfully converted to .sdf\n')
		# Delete .log files after conversion to .sdf
		[os.remove(os.path.join(babel_o[:-5], x)) for x in os.listdir(babel_o[:-5]) if x.lower().endswith('.log')]

	# Workflow for handing ORCA files: Extract ESP data, extract geometry from .out file and write a .xyz, then convert the .xyz to .sdf. Babel doesn't like ORCA .out files. 
	if orca_logs:
		print('Extracting ESP info from ORCA .out files and converting the .out files to .xyz files for processing via OpenBabel.\n')

		for file in logs:
			with open(os.path.join(directory, file), 'r') as opf:
				data_XYZ = opf.read()
				data_ESP = data_XYZ.replace('\n', 'ggez')

			#Extract ESP charges and append to ESP list
			try:
				ESP_charges = re.findall(r'CHELPG Charges(.*?)Total charge:', data_ESP)[-1]
				ESP_charges = data_ESP.split('ggez')
				ESP_charges = [x.split()[-1] for x in data_ESP if ':' in x]
				ESP.append(ESP_charges)
			except:
				error_popup('critical','ESP Charge Error',f'{file} is missing ESP data, did it finish correctly?\nGUI will now exit.')
				return

			# Create a .xyz file from ORCA .out files
			geometry = []
			GEOM = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data_XYZ)

			#parse geometry info from output
			for i in GEOM[-1].split('\n')[2:-3]:
				i = i.split()
				if len(i) > 0:
					geometry.append(i)

			#write geometry info to a formatted string
			fs = ''
			for i in geometry:
				fs += f'{i[0]:2s}    {i[1]:12s}    {i[2]:12s}    {i[3]:12s}\n'

			#write formatted geom string to a .xyz file with a header composed of the number of atoms\n filename\n xyz coordinates
			with open(os.path.join(directory, file[:-4] + '.xyz'), 'w') as opf:
				opf.write(f'{len(geometry)}\n')
				opf.write(f'{file[:-4]}\n')
				opf.write(fs + '\n')

		# Convert .xyz fto .sdf via babel (Babel has trouble converting ORCA .out files to .sdf, so we need to use this workaround)
		print('Converting ORCA xyz to sdf.\n')
		babel_i = os.path.join(directory, '*.xyz')
		babel_o = os.path.join(directory, '*.sdf')
		command = f'babel "{babel_i}" -osdf "{babel_o}"'

		try:
			convert_babel = subprocess.check_output(command, shell=True)
		except subprocess.CalledProcessError:
			try:
				print(f'Subproccess call to OpenBabel failed. Attempting to execute {command} via os.system')
				os.system(command) #alternative to subproccess call in case it fails.		
			except:
				error_popup('critical','OpenBabel Conversion Error','Encountered an error using Babel to convert the ORCA .out files to .sdf files. Is Babel installed?')
				return

		# Delete .xyz and .out files after conversion to .sdf
		[os.remove(os.path.join(babel_o[:-5], x)) for x in os.listdir(babel_o[:-5]) if x.lower().endswith(('.xyz', '.out'))]

	#Convert .sdf files to MobCal-MPI inputs
	sdf_files = [x for x in os.listdir(babel_o[:-5]) if x.endswith('.sdf')]

	for file in sdf_files:
		with open(os.path.join(babel_o[:-5], file), 'r') as opf:
			data = opf.readlines()

		#check for blackslash in first line of .sdf file and remove if present, as this will mess up file naming later on
		if data[0].find('\\') != -1:
			data[0] = data[0][data[0].rfind('\\') + 1:-5] + '\n'

		with open(os.path.join(babel_o[:-5], file), 'w') as opf:
			for line in data:
				opf.write(line)

		#sdf2tinkerxyz
		command = f'"{sdf2xyz2sdf_Directory}" < {file}'
		try:
			subprocess.check_output(command, cwd=babel_o[:-5], shell=True)
		except subprocess.CalledProcessError:
			try:
				print(f'Subproccess call to sdf2tinkerxyz failed. Attempting to execute {command} via os.system')
				os.system(command) #alternative to subproccess call in case it fails.

			except:
				error_popup('critical','sdf2tinkerxyz Error','Encountered an error using sdf2tinkerxyz.exe. Please check that the sdf2xyz2sdf directory dialog box entry contains the sdf2tinkerxyz.exe is installed there and that the sdf2xyz2sdf directory is added to your system\'s PATH (see manual).')
				return

		time.sleep(0.4) # short delay because sdf2tinkerxyz conversion is not instantaneous. 

		key_file = os.path.join(babel_o[:-5], data[0][:-1] + '.key')
		xyz_file = os.path.join(babel_o[:-5], data[0][:-1] + '.xyz')
		mfj_file = os.path.join(babel_o[:-5], data[0][:-1] + '.mfj')

		xyz_to_mfj(os.getcwd() + '\\mfj_creator\\Python\\', xyz_file, key_file, mfj_file, charge, parameters)

		#After .mfj file is created, remove the .sdf file, .key file, and tinker xyz file as they are no longer useful. 
		os.remove(os.path.join(babel_o[:-5], file))
		os.remove(key_file)
		os.remove(xyz_file)

	print(f'Process completed in {round((time.time() - start_time) / 60, 2)} minutes.')
