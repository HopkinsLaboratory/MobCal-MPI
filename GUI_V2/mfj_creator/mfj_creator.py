import os
from mfj_creator.Python.Main import run
def mfjc(self):
    errors = 0
    #Check if log list provided
    directory = self.ui.t1le_1.text().strip()
    if os.path.isdir(directory) == False:
        self.ui.t1le_1.setStyleSheet("background-color: rgb(255, 0, 0)")
        errors+=1
    #Check if log list provided
    csv_list = self.ui.t1le_2.text().strip()
    #If it isnt empty
    if len(csv_list) != 0:
        #Check if user specified path to file
        if os.path.isfile(csv_list) == False:
            #Check if user specified filename in directory
            csv_list = os.path.join(directory,csv_list)
            if os.path.isfile(csv_list) == False:
                self.ui.t1le_2.setStyleSheet("background-color: rgb(255, 0, 0)")
                errors+=1                
    #Check that the exe exists
    sdf2xyz2sdf_Directory = self.ui.t1le_3.text().strip()
    if os.path.isfile(sdf2xyz2sdf_Directory) == False:
        self.ui.t1le_3.setStyleSheet("background-color: rgb(255, 0, 0)")
        errors+=1
    #Check the temperature field for valid input
    ttemps = self.ui.t1le_4.text().split(',')
    #Check if single index is a valid integer
    if len(ttemps) == 1: # only T_bath given
        try: 
            temps = str(int(ttemps[0]))
        except ValueError: 
            self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 0, 0)")
            errors+=1
    #Check if all indexes is are valid integers
    elif len(ttemps) == 2 or len(ttemps) == 3: # T_bath and Teff_max [+ NTgrid] given; added in V2
        try: temps = ' '.join([str(int(x)) for x in ttemps])
        except ValueError:
            self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 0, 0)")
            errors+=1
        #Check if 2nd temp is larger than 1st one
        if int(ttemps[1]) < int(ttemps[0]):
            self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 0, 0)")
            errors+=1
            raise ValueError('Teff_max needs to be larger than T_bath!')
        #Check if 3rd entry is non-negative
        if len(ttemps) == 3 and int(ttemps[2]) < 0:
            self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 0, 0)")
            errors+=1
            raise ValueError('NTgrid cannot be negative!')
    else: # more than 3 entries is not allowed anymore
        self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 0, 0)")
        errors+=1
        raise ValueError('Maximum three values allowed!')
    #Check if the number of cores is correct
    v_int = self.ui.t1sb2.value()
    b_int = self.ui.t1sb3.value()
    n_cores = self.ui.t1sb4.value()
    if (v_int/n_cores).is_integer() != True:
        self.ui.t1sb4.setStyleSheet("background-color: rgb(255, 0, 0)")
        errors+=1
    if (b_int/n_cores).is_integer() != True:
        self.ui.t1sb4.setStyleSheet("background-color: rgb(255, 0, 0)")
        errors+=1
    if errors != 0 : return
    #parameters = [self.ui.t1sb1.value(),v_int,b_int,str(self.ui.t1cb2.currentIndex()+1),temps] # V1
    parameters = [self.ui.t1sb1.value(),v_int,b_int,str(self.ui.t1cb2.currentIndex()+1),
                  str(self.ui.t1cb3.currentIndex()),temps] # V1
    run(directory,csv_list,sdf2xyz2sdf_Directory,str(self.ui.t1cb1.currentText()),parameters)





