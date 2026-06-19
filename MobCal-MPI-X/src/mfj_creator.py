import os
import re
import subprocess
import time
import sys
from pathlib import Path
from shutil import copyfile
from GUI.mfj_creator.Python.xyz_to_mfj import *

def run(directory, csv, sdf2xyz2sdf_Directory, charge, parameters):
        start_time = time.time()
        sdf_files = []
        def get_files(file):
                return file.lower().endswith(('.log', '.out')) and '_atom' not in file.lower() #_atom needed for ORCA jobs that use pseudopotentials add _atom(atomic number) to the filename, which are not useful for .mfj conversion

        def estimated_runtime(number_of_files):
                #about 1.5 seconds per file, printed out in converted to minutes
                return round((number_of_files * 1.5) / (60), 2)
        
        logs = []

        if csv != '':
                with open(csv, 'r') as opf:
                        logs = [x for x in opf.read().split('\n') if get_files(x)]
        else:
                logs = [x for x in os.listdir(directory) if get_files(x)]

        print(f'Process estimated to take {estimated_runtime(len(logs))} minutes.');

        directory = directory.rstrip('/') + '/'

        #make MobCal_inputs directory (unless it already exists)
        try:
                os.mkdir(os.path.join(directory,'Mobcal_Inputs'))
        except FileExistsError:
                #if directory already exists, delete all existing contents
                [os.remove(os.path.join(os.path.join(directory, 'Mobcal_Inputs'), file)) for file in os.listdir(os.path.join(directory, 'Mobcal_Inputs'))]
        
        # Check if both Gaussian and ORCA files are present.
        gaussian_logs = any(log.lower().endswith('.log') for log in logs)
        orca_logs = any(log.lower().endswith('.out') for log in logs)

        if gaussian_logs and orca_logs:
                print('critical','File Error','Gaussian files (.log) and ORCA files (.out) are both present in the specified directory. Please put the Gaussian and ORCA files into separate directories and re-run the mfj creator.')
                return
        elif not gaussian_logs and not orca_logs:
                print('critical','File Error','There are no Gaussian files (.log) or ORCA files (.out) in the specified directory')
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
                        print('critical','OpenBabel Error',f'Open Babel version: {openbabel_version} was found on your PC. Please uninstall it, and install v2.4.1 as recommended in the manual.')
                        return
                
        except subprocess.CalledProcessError as e:
                print(f'Subprocess error when calling babel: {e.output}')
                return

        except FileNotFoundError:
                print('critical','OpenBabel Error','Open Babel could not be found on your PC. Please install v2.4.1 as stated in the manual and ensure the babel.exe is available in your system\'s PATH')
                return

        #check for missing python files
        required_python_files = os.path.dirname(os.path.join(os.path.dirname(__file__),"GUI", "mfj_creator", "Python", "test.py"))
        required_files = ['xyz_to_mfj.py', 'mass.prm', 'vdW.prm']
        missing_files = []

        for file in required_files:
                file_path = os.path.join(required_python_files, file)
                print(file_path)
                if not os.path.isfile(file_path): 
                        missing_files.append(file)

        if missing_files: #checks if list is populated. If populated, will return true. If empty, will return false. 
                missing_files_text = ', '.join(missing_files)
                print('critical','File Dependency Error',f'The following required Python files are missing from {required_python_files} : {missing_files_text}. Please redownload from the GitHub Repo and do not remove anything.')
                return

        print(f'If any errors are encountered, they will be written to: {directory}Errors.csv\n')
        print(f'To resolve the error message "No mass and vdw info for atom label XXX," please follow these steps:')
        print(f'        1. Open the .mfj file that triggered the error in a text editor.')
        print(f'        2. Locate the missing mass and van der Waals (vdW) information for the atom labeled as "XXX."')
        print(f'        3. You can find the required data in the mass.prm and vdw.prm files, which are located in the /mfj_creator/Python directory.')
        print(f'        4. Look for a similiar atom type, and fill in the missing mass and vdW values with these new ones in your .mfj file.\n')
        
        if os.path.isfile(os.path.join(directory, 'Errors.csv')): 
                os.remove(os.path.join(directory, 'Errors.csv')) #remove errors.csv if it already exists

        ESP = [] #empty list for storing ESP charges from each output file
        
        # Workflow for handing ORCA files: Extract ESP data, extract geometry from .out file and write a .xyz, then convert the .xyz to .sdf. Babel doesn't like ORCA .out files. 
        if orca_logs:
                print('Extracting ESP info from ORCA .out files and converting the .out files to .xyz files for processing via OpenBabel.\n')

                for file in logs:
                        print("Processing:", file)
                        with open(os.path.join(directory, file), 'r') as opf:
                                data_XYZ = opf.read()
                                data_ESP = data_XYZ.replace('\n', 'ggez')

                        #Extract ESP charges and append to ESP list (ESP data was not being extracted properly in previous version of main.py. Data was being processed on the native data_ESP variable, not ESP_charges. Sorry users!)
                        try:
                                ESP_charges = re.findall(r'CHELPG Charges(.*?)Total charge:', data_ESP)[-1]
                                ESP_charges = [x.split()[-1] for x in ESP_charges.split('ggez') if ':' in x]
                                ESP.append(ESP_charges)
                        except:
                                print('critical','ESP Charge Error',f'{file} is missing ESP data, did it finish correctly? MFJ file will not be generated.')
                                ESP.append([])

                        # Create a .xyz file from ORCA .out files
                        geometry = []
                        GEOM = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data_XYZ)

                        if(len(GEOM) > 0):
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
                                temp = os.path.join(directory, file[:-4] + '.xyz')
                                with open(temp, 'w') as opf:
                                        opf.write(f'{len(geometry)}\n')
                                        opf.write(f'{file[:-4]}\n')
                                        opf.write(fs + '\n')

                                # Convert .xyz fto .sdf via babel (Babel has trouble converting ORCA .out files to .sdf, so we need to use this workaround)
                                #print('Converting ORCA xyz to sdf.\n')

                                baseFileName = file.split('.')[0]
                                babel_i = os.path.join(directory, baseFileName + '.xyz')
                                babel_o = os.path.join(directory, baseFileName + '.sdf')
                                command = f'babel "{babel_i}" -osdf "{babel_o}"'

                                try:
                                        convert_babel = subprocess.check_output(command, shell=True)
                                        #Build a list of all the file names here as os.listdir returns a different order for some reason...
                                        sdf_files = sdf_files + [baseFileName + '.sdf']
                                except subprocess.CalledProcessError:
                                        try:
                                                print(f'Subproccess call to OpenBabel failed. Attempting to execute {command} via os.system')
                                                os.system(command) #alternative to subproccess call in case it fails.                
                                        except:
                                                print('critical','OpenBabel Conversion Error','Encountered an error using Babel to convert the ORCA .out files to .sdf files. Is Babel installed?')
                                                return

                # Delete .xyz and .out files after conversion to .sdf
                [os.remove(os.path.join(directory, x)) for x in os.listdir(directory) if x.lower().endswith(('.xyz', '.out'))]
                

        #Variable to track how many SDF files had errors
        errorCount = 0
        
        #Convert .sdf files to MobCal-MPI inputs        
        for file_num,file in enumerate(sdf_files):
                with open(os.path.join(directory, file), 'r') as opf:
                        data = opf.readlines()

                #check for blackslash in first line of .sdf file and remove if present, as this will mess up file naming later on
                if data[0].find('\\') != -1:
                        data[0] = data[0][data[0].rfind('\\') + 1:-5] + '\n'

                with open(os.path.join(directory, file), 'w') as opf:
                        for line in data:
                                opf.write(line)
                                
                #sdf2tinkerxyz
                command = f'"{sdf2xyz2sdf_Directory}" < {file}'
           
                try:
                        process = subprocess.check_output(command, cwd=directory, shell=True, stderr=subprocess.PIPE)
                except subprocess.CalledProcessError:
                        try:
                                print(f'Subproccess call to sdf2tinkerxyz failed. Attempting to execute {command} via os.system')
                                os.system(command) #alternative to subproccess call in case it fails.

                        except:
                                print('critical','sdf2tinkerxyz Error','Encountered an error using sdf2tinkerxyz.exe. Please check that the sdf2xyz2sdf directory dialog box entry contains the sdf2tinkerxyz.exe is installed there and that the sdf2xyz2sdf directory is added to your system\'s PATH (see manual).')
                                return

                if (len(ESP[file_num]) > 0 ):
                        #Update data in the .key file with the ESP data that was extracted earlier
                        key_filename = os.path.join(directory, data[0][:-1] + '.key')                        
                        with open(os.path.join(directory,data[0][:-1]+'.key'),'r') as key:
                                key_data = key.readlines()
                        try:

                                with open(key_filename, 'w') as key:
                                        for index, line in enumerate(key_data): # Iterate through the lines in key_data
                                                if line.split()[0] == 'charge': 
                                                        # If the line starts with charge, update the third element with ESP data
                                                        line_parts = line.split()
                                                        line_parts[2] = ESP[file_num][index]
                                                        key.write('%s %s %s\n' % (line_parts[0], line_parts[1], line_parts[2]))

                                                else:
                                                        key.write(line) # If the line does not start with charge, write the line as it is

                        except PermissionError:
                                print('critical','sdf2tinkerxyz Error',f'Cannot access: {data[0][:-1]}.key. Please restart the program.')
                                return

                        xyz_file = os.path.join(directory, data[0][:-1] + '.xyz')
                        mfj_file = os.path.join(directory, data[0][:-1] + '.mfj')

                        #xyz_to_mfj(os.getcwd() + '\\mfj_creator\\Python\\', xyz_file, key_file, mfj_file, charge, parameters)
                        xyz_to_mfj(required_python_files,directory, xyz_file, key_filename, mfj_file, charge, parameters)

                        #After .mfj file is created, remove the .sdf file, .key file, and tinker xyz file as they are no longer useful. 
                        os.remove(os.path.join(directory, file))
                        os.remove(key_filename)
                        os.remove(xyz_file)
                else:
                        errorCount += 1
                        
        print(f'Process completed in {round((time.time() - start_time) / 60, 2)} minutes.')
        print(f'Success Rate: {round((1 - errorCount/len(sdf_files))*100)}')
        
#If provided all five arguments
if (len(sys.argv) > 4):

        csvDirectory = sys.argv[2]
        if(csvDirectory == '.'):
                csvDirectory = ''

        paramList = sys.argv[4].split(',')

        #Run program
        run(sys.argv[1], "", sys.argv[2], sys.argv[3], paramList)
else:
        print("Directory, csv, sdf2xyz2sdf_Directory, charge, parameters are required inputs...")
        sys.exit(1)
        
