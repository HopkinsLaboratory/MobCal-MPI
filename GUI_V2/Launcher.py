'''
Accompanying GUI for MobCal-MPI (v2.0)
Last update: 2024-01-19

Authors:
GUI Framework: CI and AH
mjf creator: CI
mout Analyzer: AH
'''

#Import modules needed for preparation of GUI launch
import importlib
from pathlib import Path
import platform
import os
import random
import sys
import shutil 
import subprocess
import stat
import time

# Before the GUI launches, check that the user has the required packages to run the MobCal-MPI GUI
#The most troublesome package is Git, which also requires GitHub desktop to be on the user's machine. First, we check if it is installed.

def find_github_desktop():
    '''A function to dynamically locate GitHub Desktop'''
    os_type = platform.system()
    paths_to_check = []

    #Check process for Windows users
    if os_type == 'Windows':
        
        # First, try to find the path in the registry
        try:
            #Irrespective of install location of Git desktop, there should always be a registry key here for windows users
            import winreg
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r'Software\Microsoft\Windows\CurrentVersion\Uninstall\GitHubDesktop') as key:
                path, _ = winreg.QueryValueEx(key, 'InstallLocation')
                git_exe = Path(path) / 'git.exe'
                if git_exe.exists():
                    return str(git_exe.parent)
        except FileNotFoundError:
            pass
        
        #If unsucessful, try some other default locations.  
        paths_to_check.extend([
            Path(os.environ.get('LOCALAPPDATA', '')) / 'GitHubDesktop',
            Path(os.environ.get('PROGRAMFILES', '')) / 'GitHub Desktop',
            Path(os.environ.get('PROGRAMFILES(X86)', '')) / 'GitHub Desktop',
            Path(os.environ.get('USERPROFILE', '')) / 'AppData' / 'Local' / 'GitHubDesktop'
            # Users can add more Windows-specific paths here if they installed GitHub Desktop to a non-default location.
        ])

    # Check process for Mac users (I don't have a MAc so I have not been able to check if this works; I'm going off of stack exchange here)
    elif os_type == 'Darwin':
        paths_to_check.extend([
            Path('/Applications/GitHub Desktop.app'),
            Path.home() / 'Applications' / 'GitHub Desktop.app'
        ])
        # Users can add more Mac-specific paths here if they installed GitHub Desktop to a non-default location.

    # Check process for Linux users (I don't have have a Linux machine so I have not been able to check if this works; I'm going off of stack exchange here)
    elif os_type == 'Linux':
        paths_to_check.extend([
            Path('/usr/bin/github-desktop'),
            Path('/usr/local/bin/github-desktop')
        ])
        # Users can add more Mac-specific paths here if they installed GitHub Desktop to a non-default location.
    
    for path in paths_to_check:
        if (os_type == 'Windows' and path.is_dir()) or (os_type in ['Darwin', 'Linux'] and path.exists()):
            return str(path)
    
    return None

def add_to_path(new_path):
    '''Takes a path as input and adds it to the system's PATH environment variable'''
    os.environ['PATH'] = new_path + ';' + os.environ['PATH']

def check_github_desktop():
    '''Checks if GitHub Desktop is installed on the Users PC and adds it to PATH, and if it isn't, promts them to instal it before continuing'''
    
    git_desktop_path = find_github_desktop()
    if not git_desktop_path:
        print('GitHub Desktop not detected. Please install from the following URL:')
        print('https://desktop.github.com/')
        sys.exit(0)
    else:
        #If GitHub desktop is found, add it to the system's PATH
        add_to_path(git_desktop_path)

#Now with those functions set up, we can check for the required python modules 

def check_python_packages(required_packages):
    for package in required_packages:
        try:
            importlib.import_module(package)
        except ModuleNotFoundError:
            offer_package_install(package)

def offer_package_install(package):
    print(f'{package} cannot be found. Would you like to install it (y/n)?')
    if input().lower() == 'y':
        print(f'Installing {package}...')
        #git imports as gitpython despite is being called git (github why????????????????????)
        package_name = 'gitpython' if package == 'git' else package
        subprocess.run([sys.executable, '-m', 'pip', 'install', package_name])
    else:
        sys.exit('Exiting: Required package not installed.')

if __name__ == '__main__':
    #Check that github desktop is installed
    check_github_desktop()

    #check that all required python modules are installed
    required_packages = ['PyQt6', 'numpy', 'scipy', 'matplotlib', 'git']
    check_python_packages(required_packages)

#Now that all the required packages are installed, we can import the modules/functions used by the GUI
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6 import QtGui
from PyQt6.QtWidgets import QMessageBox

#import numpy as np 
#import scipy, matplotlib, git

import git
import numpy as np
from math import factorial
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib

#import GUI functions
from gui.Mobcal import Ui_Dialog
from gui.Update import Update_GUI_files
from mfj_creator.mfj_creator import mfjc
from output_analyzer.mout_Analyser import mout_info, many_mout

#The GUI is likely to the updated throughout the years, so its best practice to implement some update functionality - users may not check GitHub frequently. 
def get_latest_commit_sha(repo_url, branch='HEAD'):
    '''Grabs the SHA value associated with the latest commit to a GitHub repo. Function takes a URL as input.''' 
    try:
        # Run the git ls-remote command
        result = subprocess.run(['git', 'ls-remote', repo_url, branch], capture_output=True, text=True)

        # Check if the command was successful
        if result.returncode != 0:
            raise Exception(f'Error: {result.stderr}')

        # Parse the output to get the SHA
        output = result.stdout.split()
        return output[0] if output else None
    
    except Exception as e:
        raise Exception(f'An error occurred: {e}')

def delete_dir(directory):
    '''Windows has special permissions on read-only folders - here's a function that deals with it using shutil.rmtree'''

    def handleRemoveReadonly(func, path, exc_info):
        '''A function to deal with the error encourntered when trying to delete a file using shutil.rmtree that is read-only in Windows.'''
        
        #from: https://stackoverflow.com/questions/4829043/how-to-remove-read-only-attrib-directory-with-python-in-windows
        os.chmod(path, stat.S_IWRITE)
        os.unlink(path)

    try:
        shutil.rmtree(directory, onexc=handleRemoveReadonly)
    
    except PermissionError:
        print(f'Encountered permission error when attempting to delete {directory}. This is probably caused by a OneDrive sync issue - you can manually delete the temp folder after the GUI loads. ')
        pass

    except Exception as e: 
        raise Exception(f'Error encountered trying to remove {directory}: {e}') 

class AppWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('cleanlooks'))
       
        ## Tab 1
        self.ui.t1_pB_submit.clicked.connect(lambda: self.submit_button())
        
        ## Tab 3
        self.ui.t3_pB_load.clicked.connect(lambda: self.load_single_button()) # V2
        self.ui.t3_pB_CCSplot.clicked.connect(lambda: self.CCSplot_button()) # V2
        self.ui.t3_pB_Qlplot.clicked.connect(lambda: self.Qlplot_button()) # V2
        self.ui.t3_pB_alphaplot.clicked.connect(lambda: self.alphaplot_button()) # V2
        self.ui.t3_pB_CCSexport.clicked.connect(lambda: self.CCSexport_button()) # V2
        self.ui.t3_pB_Qlexport.clicked.connect(lambda: self.Qlexport_button()) # V2
        self.ui.t3_pB_alphaexport.clicked.connect(lambda: self.alphaexport_button()) # V2
        self.ui.t3_pB_summaryexport.clicked.connect(lambda: self.summaryexport_button()) # V2
        
        ## Tab 4
        self.ui.t4_pB_load.clicked.connect(lambda: self.load_many_button()) # V2
        self.ui.t4_pB_CCSexport.clicked.connect(lambda: self.many_CCSexport_button()) # V2
        self.ui.t4_pB_Kexport.clicked.connect(lambda: self.many_Kexport_button()) # V2
        
        ## stuff
        self.M = 0 # fake class instance later to be replaced
        self.MMs = 0 # fake class instance later to be replaced
        self.show()

    def error_popup(self, urgency, type, message):
        '''Sends a popup to the GUI if proc'd. Usage is urgency level, error title, and message text'''
        msg = QMessageBox()
        msg.setWindowTitle(type)
        msg.setText(message)

        if urgency == 'warning' or urgency == 'Warning':
            msg.setIcon(QMessageBox.Icon.Warning)
        else:
            msg.setIcon(QMessageBox.Icon.Critical)

        msg.exec()

    def reset_fields(self):
        self.ui.t1sb4.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t1le_1.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t1le_2.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t1le_3.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t1le_4.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t3le_1.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t3le_2.setStyleSheet('background-color: rgb(255, 255, 255)')
        self.ui.t4le_1.setStyleSheet('background-color: rgb(255, 255, 255)')

    ## Tab 1
    def submit_button(self):
        self.reset_fields()
        mfjc(self)

    ## Tab 3
    def load_single_button(self):
        '''Load mout file button'''
        self.reset_fields()
        directory = self.ui.t3le_1.text().strip()
        file = self.ui.t3le_2.text().strip()
        files = os.listdir(directory)
        if os.path.isdir(directory) == False:
            self.ui.t3le_1.setStyleSheet('background-color: rgb(255, 0, 0)')
            self.error_popup('critical','File path error','The .mout Directory does not exist. Please check that you have specified the correct directory.')
            return
        if not any(filename.endswith('.mout') for filename in os.listdir(directory)):
            self.ui.t3le_1.setStyleSheet('background-color: rgb(255, 0, 0)')
            self.error_popup('critical','File dependency error','The .mout directory does not contain any .mout files!')
            return

        if directory[-1] != r'/':
            directory += '/' #ensure directory ends with a /

        self.M = mout_info(directory, file) # load in class instance

        # activate buttons
        self.ui.t3_pB_CCSplot.setEnabled(True)
        self.ui.t3_pB_Qlplot.setEnabled(True)
        self.ui.t3_pB_CCSexport.setEnabled(True)
        self.ui.t3_pB_Qlexport.setEnabled(True)
        self.ui.t3_pB_alphaexport.setEnabled(True)
        self.ui.t3_pB_summaryexport.setEnabled(True)

        # if no high-field, deactivate the 'plot Mobility', export CCS, and export mobility buttons
        if self.M.highfield:
            self.ui.t3_pB_alphaplot.setEnabled(True)
        else:
            self.ui.t3_pB_alphaplot.setEnabled(False)
            self.ui.t3_pB_CCSexport.setEnabled(False)
            self.ui.t3_pB_alphaexport.setEnabled(False)

        # print short summary in the text field
        self.ui.t3le_3.setText('\n'.join(self.M.summary_text))

    def CCSplot_button(self):
        '''plots the CCS integrad (low field) or CCS-v-Teff (high field)'''
        if self.M.highfield:
            self.M.plot_CCS()
        else:
            self.M.plot_CCS_integrand()

    def alphaplot_button(self):
        '''plots the alpha function and adds the coefficients to the Summary text field'''
        if self.M.highfield:
            alpha = self.M.get_alpha_coeff(nord=6, plotting=True)
            sum_text_update = self.M.summary_text
            sum_text_update.append('')
            sum_text_update.append('alpha coefficients:')
            sum_text_update.append('a2 = %.2e Td^-2' % (alpha[0]))
            sum_text_update.append('a4 = %.2e Td^-4' % (alpha[1]))
            sum_text_update.append('a6 = %.2e Td^-6' % (alpha[2]))
            self.ui.t3le_3.setText('\n'.join(sum_text_update))
        else:
            print('no high field data')

    def Qlplot_button(self):
        '''plots the momentum transfer cross sections (l=1,2,3) over g'''
        try:
            self.M.plot_Qldat()
        except:
            self.error_popup('critical','File dependency error','No .mout file has been loaded in yet and/or the loading process has not finished!')
            return

    def CCSexport_button(self):
        '''Writes the CCS data to a file'''
        self.M.export_CCS()

    def Qlexport_button(self):
        '''Writes the CCS data to a file'''
        self.M.export_Ql()

    def alphaexport_button(self):
        '''Writes the CCS data to a file'''
        self.M.export_K()

    def summaryexport_button(self):
        '''Writes the CCS data to a file'''
        self.M.export_summary()

    ## Tab 4
    def load_many_button(self):
        '''Load mout file button'''
        self.reset_fields()
        directory = self.ui.t4le_1.text().strip()
        if os.path.isdir(directory) == False:
            self.ui.t4le_1.setStyleSheet('background-color: rgb(255, 0, 0)')
            self.error_popup('critical','File path error','The .mout Directory does not exist. Please check that you have specified the correct directory.')
            return
        if not any(filename.endswith('.mout') for filename in os.listdir(directory)):
            self.ui.t4le_1.setStyleSheet('background-color: rgb(255, 0, 0)')
            self.error_popup('critical','File dependency error','The .mout directory does not contain any .mout files!')
            return

        if directory[-1] != r'/':
            directory += '/' #ensure directory ends with a /

        self.MMs = many_mout(directory) # load in class instance
        # activate buttons
        self.ui.t4_pB_CCSexport.setEnabled(True)
        self.ui.t4_pB_Kexport.setEnabled(True)
        if self.MMs.dat_type == 'highfield':
            self.ui.t4_cb_export.setEnabled(True)
        else:
            self.ui.t4_cb_export.setCurrentIndex(0) # set back to low-field
            self.ui.t4_cb_export.setEnabled(False)

        # print a short summary in the text field
        self.ui.t4le_3.setText('\n'.join(self.MMs.summary_text))

    def many_CCSexport_button(self):
        '''Writes the CCS data to a file'''
        itype = self.ui.t4_cb_export.currentIndex() # low or high field export
        self.MMs.export_CCS(itype)

    def many_Kexport_button(self):
        '''Writes the CCS data to a file'''
        itype = self.ui.t4_cb_export.currentIndex() # low or high field export
        self.MMs.export_Mobility(itype)

    def check_for_update(self):
        
        # Get the current working directory and define the temporary directory path
        root = os.getcwd()

        #URL of the MobCal-MPI repo
        repo_url = 'https://github.com/ChristianIeritano/MobCal-MPI-UpdateTesting'

        #get SHA value of repo-ID
        repo_SHA = get_latest_commit_sha(repo_url)

        #get SHA value of local repo

        # File to store the local version's SHA value
        ID_file = os.path.join(root, 'ID.txt') 
        try:
            with open(ID_file,'r') as opf:
                file_content = opf.read()

                # Check if the file is in the correct format
                if '\n' in file_content or '\r' in file_content or ' ' in file_content:
                    raise Exception(f'{ID_file} is not in the correct format. Please re-download this file from the GitHub repo and re-run the GUI launcher.')

                local_SHA = file_content.strip()

        except FileNotFoundError:
            raise FileNotFoundError(f'{ID_file} could not be found. Please re-download from the MobCal-MPI GitHub repo and re-run the GUI launcher.')
            
        # Remove the temporary directory if it exists and only if the local and repo SHAs match
        if repo_SHA == local_SHA: 
            temp_dir = os.path.join(os.getcwd(), 'temp')
            if os.path.isdir(temp_dir):
                delete_dir(temp_dir)

        else:
            choice_title = 'Update available!'
            choice_prompt = 'An update to the MobCal-MPI is available. Note that this can include updates to the MobCal-MPI source code and/or the GUI, so please check Github for details.\nWould you like to update now?'
            choice = QMessageBox.question(self, choice_title, choice_prompt, QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)

            if choice == QMessageBox.StandardButton.Yes:
                print('The MobCal-MPI GUI is being updated. Any errors encountered during the update process will be printed below.')
                
                #Run the update function
                Update_GUI_files(repo_url, root, ID_file, repo_SHA, delete_dir)

                print('The MobCal-MPI GUI files have been succesfully updated to their current version. The GUI will now close. Please reload the GUI')
                sys.exit(0)

            else:
                print('The user has opted to use their local version of MobCal-MPI.')

    def close_application(self):
        choice_title = 'Exit Confirmation'
        choice_prompt = 'Are you sure you wish to exit?'
        choice = QMessageBox.question(self, choice_title, choice_prompt, QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)

        if choice == QMessageBox.StandardButton.Yes:
            sys.exit(0)

    def closeEvent(self, event):
        event.ignore() # Do not let the application close without a prompt
        self.close_application()

def main():
    app = QtWidgets.QApplication(sys.argv)
    main_window = AppWindow()
    main_window.show()  # Show the main window
    main_window.check_for_update()  # Check for updates after showing the window        
    sys.exit(app.exec())  # Changed from exec_() to exec() - PyQT6

sys._excepthook = sys.excepthook

def exception_hook(exctype, value, traceback):
    print(exctype, value, traceback)
    sys._excepthook(exctype, value, traceback)

sys.excepthook = exception_hook

if __name__ == '__main__':
    main()

