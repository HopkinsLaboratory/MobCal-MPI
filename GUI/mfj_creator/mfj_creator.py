import os
from mfj_creator.Python.Main import run
from PyQt6.QtWidgets import QMessageBox

def mfjc(self):

	# Check if a valid directory is provided
	directory = self.ui.t1le_1.text().strip()
	if not os.path.isdir(directory):
		self.ui.t1le_1.setStyleSheet('background-color: rgb(255, 0, 0)')
		self.error_popup('critical','File path error','The Out/Log Directory does not exist. Please check that you have specified the correct directory.')

	# Check if a valid CSV list is provided
	csv_list = self.ui.t1le_2.text().strip()
	if csv_list:
		if not os.path.isfile(csv_list):
			csv_list = os.path.join(directory, csv_list) #if directory wasn't provided for csv list, join it with the initial directory, and check if it exists again
			if not os.path.isfile(csv_list): 
				self.ui.t1le_2.setStyleSheet('background-color: rgb(255, 0, 0)')
				self.error_popup('critical','.csv file error','The Out/Log List .csv file cannot be found. Please ensure it is spelled correctly and is located within the Out/Log Directory.')

	# Check if the sdf2xyz2sdf executable exists
	sdf2xyz2sdf_Directory = self.ui.t1le_3.text().strip()
	if not os.path.isfile(sdf2xyz2sdf_Directory):
		self.ui.t1le_3.setStyleSheet('background-color: rgb(255, 0, 0)')
		self.error_popup('critical','File path error','The sdf2xyz2sdf directory specified is not valid. Please ensure that sdf2tinkerxyz.exe can be found within this directory.')

	# Check the temperature field for valid input
	ttemps = self.ui.t1le_4.text().split(',')

	# If one temperature is given, ensure that it is an integer
	if len(ttemps) == 1:
		try:
			int(ttemps[0])
		except ValueError: #if someone puts text into this block or does not enter a temperature
			self.ui.t1le_4.setStyleSheet('background-color: rgb(255, 0, 0)')
			self.error_popup('warning','Teff error','The temperature array can only be numbers and cannot be left empty.')
			return

	# If more than 3 entries are given in the temp list, ensure that all entries are integers
	elif len(ttemps) == 3:
		try:
			[int(x) for x in ttemps]
		except ValueError: # if someone puts text into this block
			self.ui.t1le_4.setStyleSheet('background-color: rgb(255, 0, 0)')
			self.error_popup('warning','Teff error','Tbath, T_eff_max, and the temperature grids must be integers separated by a comma.')
			return

		# Check if the second temperature is larger than the first one
		if int(ttemps[1]) <= int(ttemps[0]):
			self.ui.t1le_4.setStyleSheet('background-color: rgb(255, 0, 0)')
			self.error_popup('warning','Teff error','T_eff_max (entry #2) must be larger than Tbath (entry #1).')
			return

	# If the temperature list is not 1 or 3 entries long, raise an error
	else:
		self.ui.t1le_4.setStyleSheet('background-color: rgb(255, 0, 0)')
		self.error_popup('warning','Teff error','The temp list must contain either 1 entry (T_bath), or three entries (T_bath, T_eff_max, and a grid size separated by a comma). Other array sizes are not permitted. Please see the .mfj creator section of the manual for more details.')
		return

	# Check if any negative temperatures or grid spacing are given
	for ttemp in ttemps:
		if int(ttemp) <= 0:
			self.ui.t1le_4.setStyleSheet('background-color: rgb(255, 0, 0)')
			self.error_popup('warning','Teff error','Temperatures and/or grid spacing cannot be negative!')
			return

    #get values for v_int (inp), b_int (imp), and number of cores, respectively. 
	v_int = self.ui.t1sb2.value()
	b_int = self.ui.t1sb3.value()
	n_cores = self.ui.t1sb4.value()
	
	# Number of velocity grid points (inp) and impact parameters points (imp) must be divisible by the number of cores the user intends on using.
	if not (v_int / n_cores).is_integer() or not (b_int / n_cores).is_integer():
		self.ui.t1sb4.setStyleSheet('background-color: rgb(255, 0, 0)')
		self.error_popup('warning','Error','Number of velocity grid points (inp) and impact parameters points (imp) must be divisible by the number of cores that the user intends on using.')

	# If there are no errors, create the parameters list and start the process of .mfj file generation
	temps = ' '.join(map(str, ttemps))
	parameters = [
		self.ui.t1sb1.value(),
		v_int,
		b_int,
		str(self.ui.t1cb2.currentIndex() + 1),
		str(self.ui.t1cb3.currentIndex()),
		temps
	]

	# Run the 'run' function with the provided parameters
	run(self, directory, csv_list, sdf2xyz2sdf_Directory, str(self.ui.t1cb1.currentText()), parameters)
