#Specify a directory containing log files to be converted into mobcal inputs
#i.e directory = r'C:\Users\McMahon003\Desktop\Mobcal Script\Test Files'
directory = r'G:\Hopkins_Laboratory\Lipidz_Tingz\Ag_MobCalMPI_Cal\ChelpG\Outputs'

#The csv_list will read specified log files in  a .csv
#i.e csv_list = r'C:\Users\McMahon003\Desktop\Mobcal Script\Test Files\Logs.csv'
#If no files are to be specified set csv_list = r''
csv_list = r''

#Specify the directory that sdf2xyz2sdf was installed in
#i.e sdf2xyz2sdf_Directory = r'C:\open3dtools\bin\sdf2tinkerxyz.exe'
sdf2xyz2sdf_Directory = r'C:\open3dtools\bin\sdf2tinkerxyz.exe'

#This sets the charge for the mfj file, you can specify: 'calc','equal', or 'none'. Must be wrapped in quotes.
charge = 'calc'

#Specify the number of cycles of mobility calculations to run. Recommended value is 10. Must be an integer.
cycles = 10

#Specify the number of velocity integration points to take. Recommended value is 48. Must be an integer.
v_integrations = 48

#Specify the number of impact parameter integration points to take. Recommended value is 512. Must be an integer.
b_integrations = 512

#Specify the number of cores you plan to run trajectory method calcs on. Not used in input file creation, but used in error checking. Must be an integer. 
n_cores = 16

#Specify the buffer gas to use. Currently, MobCal-MPI can only handle gas = 'He' OR gas = 'N2' as an input. Must be wrapped in quotes
gas = 'N2'

#After you have edited the above fields please run in idle using F5, double click file, or using cmd prompt
#cmd prompt instructions:
#Move to the directory this file is located in, i.e: "cd C:\Users\JamFlex\ModJed\MobcalMPI_Script\mfj_creator.py"
#then type in: python Edit_Me.py

####################### ERROR HANDLING #########################
############### DO NOT MOFIDY PAST THIS POINT ##################

if gas == 'N2':
        print('Creating .mfj inputs for calculation in N2\n')
elif gas == 'He':
        print('Creating .mfj inputs for calculation in He\n')
else:
        print('Unrecognized gas defined.\n')

from Python.Main import *

a = 0
try:
    if charge != 'calc' and charge != 'equal' and charge != 'none':
        a = 1
        print('You have entered an invalid charge parameter, please see comment on line 14 for more information.')
except NameError:
    print('You have deleted the charge variable, please replace on line 15')

try:
    if os.path.isfile(sdf2xyz2sdf_Directory) != True:
        a = 1
        print('The specified sdf2xyz2sdf Directory is incorrect.')
except NameError:
    a = 1
    print('You have deleted the sdf2xyz2sdf_Directory variable, please replace on line 17')

try:
    if os.path.isdir(directory) != True:
        a = 1
        print('Directory is invalid')
except NameError:
    a = 1
    print('You have deleted the directory variable, please replace on line 3')

parameters = [cycles,v_integrations,b_integrations,gas]

if cycles and type(cycles) != int:
    a = 1
    print('n_cycles is not an integer. Please fix on line XX')

if n_cores and type(n_cores) != int:
    a = 1
    print('n_cores is not an integer. Please fix on line 30')

try:
    if v_integrations and type(v_integrations) != int:
        a = 1
        print('v_integrations is not an integer!')
    elif (v_integrations/n_cores).is_integer() != True:
        a = 1
        print('v_integrations is not divisible by n_cores.')
        
except TypeError:
    a=1
    if type(n_cores) != int:
        print('n_cores is not an integer!')
    if type(v_integrations) != int:
        print('v_integrations is not an integer!')

if a == 0:
    try:
        if b_integrations and type(b_integrations) != int:
            a = 1
            print('b_integrations is not an integer!')
        elif (b_integrations/n_cores).is_integer() != True:
            a = 1
            print('b_integrations is not divisible by n_cores.')
    except TypeError:
        a = 1
        if type(n_cores) != int:
            print('n_cores is not an integer!')
        if type(b_integrations) != int:
            print('b_integrations is not an integer!')

if gas == 'N2':
    parameters[3]=2
elif gas == 'He':
    parameters[3]=1
else:
    a = 1
    print("You have specified an invalid buffer gas. Please enter 'N2' or 'He' on line XX")

try:
    if a == 0:
        if os.path.isfile(csv_list) == True:
            run(directory,csv_list,sdf2xyz2sdf_Directory,charge,parameters)
        if os.path.isfile(csv_list) == False and csv_list != '':
            print('The specified csv_list file does not exist.')
        if os.path.isfile(csv_list) == False and csv_list == '':
            run(directory,csv_list,sdf2xyz2sdf_Directory,charge,parameters)
except ValueError:
    print('You have deleted the csv_list variable, please replace on line 8')




