import os
import sys
import re
import subprocess
import time
from shutil import copyfile
from Python.xyz_to_mfj import *

def run(directory,csv,sdf2xyz2sdf_Directory,charge,parameters):
	start_time = time.time()
	
	if csv != '':
		opf = open(csv,'r')
		logs = opf.read().split('\n')
		opf.close()
		logs = [x for x in logs if x.lower().endswith('.xyz') or x.lower().endswith('out')]
	else:
		logs = [x for x in os.listdir(directory) if x.lower().endswith('.xyz') or x.lower().endswith('.out')]
	print('Process estimated to take '+str(round((len(logs)*2.5+10)/60,2))+' minutes.')
	print(logs)
	
	if directory[-1] != '\\':
		directory = directory+'\\'
	try:
		os.mkdir(directory+'Mobcal_Inputs')
	except:
		pass
	for file in logs:
		try:
			copyfile(directory+file,directory+'Mobcal_Inputs\\'+file)
			#print('copying completed')
		except (PermissionError):
			pass
	directory = directory+'Mobcal_Inputs'

	App_Data = os.listdir(os.getenv('APPDATA'))
	Babel_Installed = [x for x in App_Data if x == 'OpenBabel-2.4.1']
	if not Babel_Installed:
		print('Open Babel could not be found on your pc, are you running the version provided?')
		sys.exit()
	
	if not os.path.isdir(os.getcwd()+'\\Python'):
		print('The required python files are missing, please reinstall')
		sys.exit()		 
	
	directory = directory.replace('//','\\') #Make sure it is in command prompt lingo
	if directory[-1] != '\\': #If the last character is not a \ add a slash
		directory = directory+'\\'
	
	print('If any errors are encountered they will be written to: '+directory+'Errors.csv')
	if os.path.isfile(directory+'\Errors.csv'): #Delete old error files before each run
		os.remove(directory+'\Errors.csv')

	##### -- ESP Data -- #####
	ESP = []
	#print(directory)
	logs = [x for x in os.listdir(directory) if x.lower().endswith('.out')]
	print(logs)
	print('Extracting ESP info from logs.\n')
	for file in logs:
		opf = open(directory+file,'r')
		data = opf.read().replace('\n','ggez')
		opf.close()

		try:
			#data = re.findall(r'Mulliken charges:(.*?)Sum of Mulliken charges',data)[-1]
			data = re.findall(r'CHELPG Charges(.*?)Total charge:',data)[-1]
			#data = re.findall(r'LOEWDIN ATOMIC CHARGES(.*?)LOEWDIN REDUCED ORBITAL CHARGES',data)[-1]
			#print(data)
			data = data.split('ggez')
			#print(data)
			data = [x.split()[-1] for x in data if ':' in x]
			#print(data)
			ESP.append(data)
			#print(ESP)
		except:
			print(file+' is missing ESP data, did it finish correctly?')


	##### -- Babel -- #####
	print('Converting logs to sdf.\n')
	babel_i = directory+'*.xyz'
	babel_o = directory+'*.sdf'

	command = str('babel "'+babel_i+'" -osdf "'+babel_o+'"') #Create command
	#print(command)
        
	try:
		convert_babel = subprocess.check_output(command, shell=True) #pass to cmd prompt
	except subprocess.CalledProcessError:
		os.system(command)
		print('Encountered error opening babel with subprocess.')

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
			print('Encountered error opening sdf2xyz with subprocess.')
		time.sleep(1)

		key = open(babel_o[:-5]+data[0][:-1]+'.key','r')
		key_data = key.readlines()
		key.close()

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
			print('Cannot access: '+data[0][:-1]+'.key'+' please restart the program.')
			sys.exit()
	 
		#print('Creating '+file[:-4]+'.mfj\n')

		key_file = babel_o[:-5]+data[0][:-1]+'.key'
		xyz_file = babel_o[:-5]+data[0][:-1]+'.xyz'
		mfj_file = babel_o[:-5]+data[0][:-1]+'.mfj'

		xyz_to_mfj(os.getcwd()+r'\Python\\',xyz_file,key_file,mfj_file,charge,parameters)

		#Delete old files
		os.remove(babel_o[:-5]+file) #sdf file
		os.remove(key_file)
		os.remove(xyz_file)

		file_num+=1
	print('Process completed in '+str(round((time.time() - start_time)/60,2))+' minutes.')

