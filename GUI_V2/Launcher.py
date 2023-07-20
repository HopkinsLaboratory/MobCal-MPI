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
from PyQt5.QtWidgets import QMessageBox

from gui.Mobcal import Ui_Dialog
from mfj_creator.mfj_creator import mfjc
from output_analyzer.mout_Analyser import mout_info, many_mout

class AppWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('cleanlooks'))
        ## Tab 1
        # self.ui.pushButton.clicked.connect(lambda: self.submit_button())
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

    def reset_fields(self):
        self.ui.t1sb4.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_1.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_2.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_3.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t1le_4.setStyleSheet("background-color: rgb(255, 255, 255)")
        #self.ui.t2le_1.setStyleSheet("background-color: rgb(255, 255, 255)")
        self.ui.t3le_1.setStyleSheet("background-color: rgb(255, 255, 255)")
    
    ## Tab 1
    def submit_button(self):
        self.reset_fields()
        mfjc(self)
        # pt_index = self.ui.tabWidget.currentIndex()
        # if pt_index == 0: mfjc(self)
        # # if pt_index == 1: ew_run(self)
        # if pt_index == 1: ccs_run(self) # output read function
    
    ## Tab 3
    def load_single_button(self):
        '''Load mout file button'''
        directory = self.ui.t3le_1.text().strip()
        file = self.ui.t3le_2.text().strip()
        if os.path.isdir(directory) == False:
            self.ui.t3le_1.setStyleSheet("background-color: rgb(255, 0, 0)")
            return
        
        if directory[-1] != r'/':
            directory += '/'
        
        self.M = mout_info(directory,file) # load in class instance
        # activate buttons
        self.ui.t3_pB_CCSplot.setEnabled(True)
        self.ui.t3_pB_Qlplot.setEnabled(True)
        self.ui.t3_pB_CCSexport.setEnabled(True)
        self.ui.t3_pB_Qlexport.setEnabled(True)
        self.ui.t3_pB_alphaexport.setEnabled(True)
        self.ui.t3_pB_summaryexport.setEnabled(True)
        # if no highfield, deactivate the "plot Mobility" button
        if self.M.highfield:
            self.ui.t3_pB_alphaplot.setEnabled(True)
        else:
            self.ui.t3_pB_alphaplot.setEnabled(False)
        
        # print short summary in text field        
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
            sum_text_update.append('a2 = %.2e Td^-2' %(alpha[0]))
            sum_text_update.append('a4 = %.2e Td^-4' %(alpha[1]))
            sum_text_update.append('a6 = %.2e Td^-6' %(alpha[2]))
            self.ui.t3le_3.setText('\n'.join( sum_text_update ))
        else:
            print('no high field data')
    
    def Qlplot_button(self):
        '''plots the momentum transfer cross sections (l=1,2,3) over g'''
        try:
            self.M.plot_Qldat()
        except:
            print('No .mout file loaded in yet or loaded file did not finish!')
    
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
        directory = self.ui.t4le_1.text().strip()
        if os.path.isdir(directory) == False:
            self.ui.t4le_1.setStyleSheet("background-color: rgb(255, 0, 0)")
            return
        
        if directory[-1] != r'/':
            directory += '/'
        
        self.MMs = many_mout(directory) # load in class instance
        # activate buttons
        self.ui.t4_pB_CCSexport.setEnabled(True)
        self.ui.t4_pB_Kexport.setEnabled(True)
        if self.MMs.dat_type == 'highfield':
            self.ui.t4_cb_export.setEnabled(True)
        else:
            self.ui.t4_cb_export.setEnabled(False)
        
        # print short summary in text field        
        self.ui.t4le_3.setText('\n'.join(self.MMs.summary_text))
    
    def many_CCSexport_button(self):
        '''Writes the CCS data to a file'''
        itype = self.ui.t4_cb_export.currentIndex() # low or high field export
        self.MMs.export_CCS(itype)
    
    def many_Kexport_button(self):
        '''Writes the CCS data to a file'''
        itype = self.ui.t4_cb_export.currentIndex() # low or high field export
        self.MMs.export_Mobility(itype)
    
    ##
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
