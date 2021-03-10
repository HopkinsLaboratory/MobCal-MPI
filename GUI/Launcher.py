import importlib
import os
import sys

#Required Packages to run the script, ask user to install if not present
rp = ['PyQt5','numpy']
installed = 0
for package in rp:
    try: globals()[package] = importlib.import_module(package)
    except (ModuleNotFoundError,ImportError) as e:
        print(f'{package} cannot be found, would you like to install it (y/n)?')
        answer = input()
        if answer == 'y':
            os.system(f'pip install {package}')
            installed+=1
    if installed > 0: sys.exit()

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui

from gui.Mobcal import Ui_Dialog
from mfj_creator.mfj_creator import mfjc
from Energy_Weighter.Energy_weighter import ew_run
from CCS_extractor.CCS_Extractor import ccs_run

class AppWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('cleanlooks'))
        self.ui.pushButton.clicked.connect(lambda: self.submit_button())
        self.ui.pushButton_2.clicked.connect(lambda: self.file_open('*.csv'))
        self.show()

    def reset_fields(self):
        self.ui.t1sb4.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_1.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_2.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_3.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t2le_1.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t3le_1.setStyleSheet("background-color: rgb(255, 255, 255)")

    def submit_button(self):
        self.reset_fields()
        pt_index = self.ui.tabWidget.currentIndex()
        if pt_index == 0: mfjc(self)
        if pt_index == 1: ew_run(self)
        if pt_index == 2: ccs_run(self)

    def file_open(self, file_type):
        name = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File','',file_type)
        if len(name) > 1: #Linux stores name in tuple
            name = name[0]
        self.ui.t2le_1.setText(name)

    def close_application(self): #Exit alert for user
        choice_title = 'Exit Confirmation'
        choice_prompt = 'Are you sure you wish to exit?'
        choice = QtWidgets.QMessageBox.question(self, choice_title, choice_prompt, QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        if choice == QtWidgets.QMessageBox.Yes:
            sys.exit()

    def closeEvent(self, event):
        event.ignore() #Do not let application close without prompt
        self.close_application() #Call the exit prompt
        
def main():
    app = QtWidgets.QApplication(sys.argv)
    _ = AppWindow()
    sys.exit(app.exec_())

sys._excepthook = sys.excepthook 
def exception_hook(exctype, value, traceback):
    print(exctype, value, traceback)
    sys._excepthook(exctype, value, traceback) 
sys.excepthook = exception_hook 

main()
