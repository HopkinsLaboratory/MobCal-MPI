import os
import re
import subprocess
import time
from shutil import copyfile
from mfj_creator.Python.xyz_to_mfj import *
from PyQt5.QtWidgets import QMessageBox

def run(directory,csv,sdf2xyz2sdf_Directory,charge,parameters):
	start_time = time.time()
	
	if csv != '':
		opf = open(csv,'r')
		logs = opf.read().split('\n')
		opf.close()
		logs = [x for x in logs if x.lower().endswith('.log') or x.lower().endswith('.out') and '_atom' not in x]
	else:
		logs = [x for x in os.listdir(directory) if x.lower().endswith('.log') or x.lower().endswith('.out') and '_atom' not in x]
	print('Process estimated to take '+str(round((len(logs)*2.5+10)/(60*2),2))+' minutes.')
	
	if directory[-1] != '\\':
		directory = directory+'\\'
	try:
		os.mkdir(directory+'Mobcal_Inputs')
	except:
		pass

	#Check if both Gaussian and ORCA files are present. If so, kill the program
	check = [end in ','.join(logs) for end in ['.log','.out']]
	if check == [1,0] or check == [0,1]:	#[1,0] is only gaussian .log and [0,1] if only ORCA .out
		pass
	elif check == [1,1]:
		msg = QMessageBox()
		msg.setWindowTitle('File Error')
		msg.setText('Gaussian files (.log) and ORCA files (.out) are both present in the specificed directory. Please put the Gaussian and ORCA files into separate directories, and re-run the mfj creator.')
		msg.setIcon(QMessageBox.Critical)
		msg.exec_() #show messagebox
		return		
	else:
		msg = QMessageBox()
		msg.setWindowTitle('File Error')
		msg.setText('There are no Gaussian files (.log) or ORCA files (.out) in the specified directory')
		msg.setIcon(QMessageBox.Critical)
		msg.exec_() #show messagebox
		return
		#print('Gaussian files (.log) and ORCA files (.out and .xyz) are present in the specificed directory. Please separate them')
		#sys.exit()

	#Copy files so that the originals don't get changed	
	for file in logs:
		try:
			copyfile(directory+file,directory+'Mobcal_Inputs\\'+file)
		except (PermissionError):
			pass
	directory = directory+'Mobcal_Inputs' #Rename global directory variable to new directory

	#Other error handling

	#Check if Babel is on PC
	App_Data = os.listdir(os.getenv('APPDATA'))
	Babel_Installed = [x for x in App_Data if x == 'OpenBabel-2.4.1']
	if not Babel_Installed:
		msg = QMessageBox()
		msg.setWindowTitle('Software Error')
		msg.setText('Open Babel could not be found on your pc, are you running the version provided?')
		msg.setIcon(QMessageBox.Critical)
		msg.exec_() #show messagebox
		return
	
	if not os.path.isdir(os.getcwd()+'\mfj_creator\\Python'):
		msg = QMessageBox()
		msg.setWindowTitle('Software Error')
		msg.setText('The required python files are missing. Please redownload from the GitHub Repo and do not remove anything.')
		msg.setIcon(QMessageBox.Critical)
		msg.exec_() #show messagebox
		return
	
	directory = directory.replace('//','\\') #Make sure it is in command prompt lingo
	if directory[-1] != '\\': #If the last character is not a \ add a slash
		directory = directory+'\\'
	
	print('If any errors are encountered they will be written to: '+directory+'Errors.csv')
	if os.path.isfile(directory+'\Errors.csv'): #Delete old error files before each run
		os.remove(directory+'\Errors.csv')

	##### -- ESP Data -- #####
	ESP = []
	logs = [x for x in os.listdir(directory) if x.lower().endswith('.log') or x.lower().endswith('.out')]
	
	#Extracting ESP data from Gaussian .log files
	if check == [1,0]:
		print('Extracting ESP info from Gaussian logs.\n')

		for file in logs:
			opf = open(directory+file,'r')
			data = opf.read().replace('\n','ggez')
			opf.close()

			#Get ESP charge data from Gaussian files
			try:
				#data = re.findall(r'Mulliken charges:(.*?)Sum of Mulliken charges',data)[-1]
				data = re.findall(r'ESP charges:(.*?)Sum of ESP charges',data)[-1]
				data = data.split('ggez')
				data = [x.split()[2] for x in data if len(x.split()) == 3]
				ESP.append(data)
			except:
				msg = QMessageBox()
				msg.setWindowTitle('ESP Charge Error')
				msg.setText(file+' is missing ESP data, did it finish correctly?')
				msg.setIcon(QMessageBox.Critical)
				msg.exec_() #show messagebox
				return
	
	#Extracting ESP data from ORCA .out files
	if check == [0,1]:
		print('Extracting ESP info from ORCA .out files.\n')

		for file in logs:
			opf = open(directory+file,'r')
			data = opf.read().replace('\n','ggez')
			opf.close()		

		#Get ESP charge data from ORCA out files
			try:
				data = re.findall(r'CHELPG Charges(.*?)Total charge:',data)[-1]
				data = data.split('ggez')
				data = [x.split()[-1] for x in data if ':' in x]
				ESP.append(data)
			except:
				msg = QMessageBox()
				msg.setWindowTitle('ESP Charge Error')
				msg.setText(file+' is missing ESP data, did it finish correctly?')
				msg.setIcon(QMessageBox.Critical)
				msg.exec_() #show messagebox
				return

	#There is a problem when converting .out files to .sdf using BABEL.... So we're going to make a .xyz file from the ORCA .out because that conversion works just fine 
	outs = [x for x in os.listdir(directory) if x.lower().endswith('.out')]
	for filename in outs:
		opf = open(directory+'\\'+filename,'r')
		data = opf.read()
		opf.close()

		#Get the geometry and write to list called geoemtry
		geometry=[] #Get geometry of atoms from file
		GEOM = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)',data)

		#Loop through each line in the geom block, split each line, then append to geometry list
		for i in GEOM[-1].split('\n')[2:-3]:
			i = i.split()
			if len(i) > 0:
				geometry.append(i)

		#Write the geometry to a formatted string
		fs = ''
		for i in geometry:
			fs+='%2s    %12s    %12s    %12s\n'%(i[0],i[1],i[2],i[3])    

		#At long last, we can write all the info to the .xyz file
		opf = open(directory+'\\'+filename[:-4]+'.xyz','w')
		opf.write(str(len(geometry))+'\n')
		opf.write(str(filename[:-4])+'\n') #Write filename without extension as the title in the xyz file (Needed for conversion to sdf)
		opf.write(fs)
		opf.write('\n')
		opf.close()

	##### -- Babel -- #####
	
	# For handling Gaussian log files
	if check == [1,0]:
	
		print('Converting logs to sdf.\n')

		babel_i = directory+'*.log'
		babel_o = directory+'*.sdf'

		command = str('babel "'+babel_i+'" -osdf "'+babel_o+'"') #Create command

		try:
			convert_babel = subprocess.check_output(command, shell=True) #pass to cmd prompt
		except subprocess.CalledProcessError:
			os.system(command)
			msg = QMessageBox()
			msg.setWindowTitle('Software Error')
			msg.setText('Encountered error opening babel with subprocess.')
			msg.setIcon(QMessageBox.Critical)
			msg.exec_() #show messagebox
			return

		#Delete .log files after conversion to .sdf
		[os.remove(babel_o[:-5]+x) for x in os.listdir(babel_o[:-5]) if x.lower().endswith('.log')]
	
	# For handling ORCA .xyz files (BABEL doesn't seem to like converting the .out file to a .xyz)
	if check == [0,1]:
	
		print('Converting ORCA xyz to sdf.\n')

		babel_i = directory+'*.xyz'
		babel_o = directory+'*.sdf'

		command = str('babel "'+babel_i+'" -osdf "'+babel_o+'"') #Create command
		
		try:
			convert_babel = subprocess.check_output(command, shell=True) #pass to cmd prompt
		except subprocess.CalledProcessError:
			os.system(command)
			msg = QMessageBox()
			msg.setWindowTitle('Software Error')
			msg.setText('Encountered error opening babel with subprocess.')
			msg.setIcon(QMessageBox.Critical)
			msg.exec_() #show messagebox
			return

		#Delete .xyz and .out files after conversion to .sdf
		[os.remove(babel_o[:-5]+x) for x in os.listdir(babel_o[:-5]) if x.lower().endswith('.xyz') or x.lower().endswith('.out')]

	files = [x for x in os.listdir(babel_o[:-5]) if x.endswith('.sdf')]
	file_num = 0 #loop through by counter rather than range(len(files))

	for file in files:
		opf = open(babel_o[:-5]+file,'r')
		data = opf.readlines()
		opf.close()

		if data[0].find('\\') != -1:
			data[0] = (data[0][data[0].rfind('\\')+1:-5]+'\n') #Fix filename

		#print('Converting '+data[0][:-1]+'.sdf to .xyz and .key files')

		opf = open(babel_o[:-5]+file,'w')
		for line in data:
			opf.write(line)
		opf.close()

		##### -- sdf2tinkerxyz -- #####
		command = '"'+sdf2xyz2sdf_Directory+'" < '+file
		try:
			subprocess.check_output(command, cwd=babel_o[:-5], shell=True) #pass to cmd prompt
		except subprocess.CalledProcessError:
			os.system(command)
			msg = QMessageBox()
			msg.setWindowTitle('Software Error')
			msg.setText('Encountered error opening sdf2xyz with subprocess.')
			msg.setIcon(QMessageBox.Critical)
			msg.exec_() #show messagebox
			return
		
		#Time delay of 0.5s; sdf2xyz conversion isn't instantaneous 
		time.sleep(0.5)

		key = open(babel_o[:-5]+data[0][:-1]+'.key','r')
		key_data = key.readlines()
		key.close()

		#Replace the default MMFF94 charges in the .key file with the ones calculated from ESP methods
		try:
			key = open(babel_o[:-5]+data[0][:-1]+'.key','w')
			for index in range(len(key_data)):
				if key_data[index].split()[0] == 'charge':
					key_data[index] = (key_data[index].split())
					key_data[index][2] = ESP[file_num][index]
					key.write('%s %s %s\n'%(key_data[index][0],key_data[index][1],key_data[index][2]))
				else:
					key.write(key_data[index])
			key.close()
		except PermissionError:
			msg = QMessageBox()
			msg.setWindowTitle('File Error')
			msg.setText('Cannot access: '+data[0][:-1]+'.key'+' please restart the program.')
			msg.setIcon(QMessageBox.Critical)
			msg.exec_() #show messagebox	
			return		
	 
		#print('Creating '+file[:-4]+'.mfj\n')

		key_file = babel_o[:-5]+data[0][:-1]+'.key'
		xyz_file = babel_o[:-5]+data[0][:-1]+'.xyz'
		mfj_file = babel_o[:-5]+data[0][:-1]+'.mfj'

		#Pass files to xyz_to_mfj python code
		xyz_to_mfj(os.getcwd()+r'\mfj_creator\\Python\\',xyz_file,key_file,mfj_file,charge,parameters)

		#Delete old files
		os.remove(babel_o[:-5]+file) #sdf file
		os.remove(key_file)
		os.remove(xyz_file)

		file_num+=1
	
	print('Process completed in '+str(round((time.time() - start_time)/60,2))+' minutes.')
