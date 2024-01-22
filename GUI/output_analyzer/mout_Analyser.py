import os
import re
import numpy as np
from math import factorial
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('QtAgg') #Use matplotlib backend that is compatible w/ PyQt6 to prevent competition w/ GUI event loop
import matplotlib.pyplot as plt

from PyQt6.QtWidgets import QMessageBox

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'lines.linewidth': 1.5})

#because of class within a class layout, its easier to call the error popup function externally rather than nest it within each new class defined in this file
def error_popup(urgency, type, message):
	'''Sends a popup to the GUI if proc'd. Usage is urgency level, error title, and message text'''
	msg = QMessageBox()
	msg.setWindowTitle(type)
	msg.setText(message)

	if urgency == 'warning' or urgency == 'Warning':
		msg.setIcon(QMessageBox.Icon.Warning)
	else:
		msg.setIcon(QMessageBox.Icon.Critical)

	msg.exec()

class mout_info:
	'''class with all info and functions from mout files'''
	def __init__(self,directory,file):
		'''reading in all data from the mout file and storing them in arrays'''
		
		try:
			opf = open(directory+file,'r')
			data = opf.read()
			opf.close()
		except: #check if mout file exists
			error_popup('critical','File Error',f'Could not find the file {file}\n Does it exist in the directory specified and/or did it finish correctly?')
			return 


		# constants without self.
		elC = 1.60217733E-19 # C
		kB = 1.380658E-23 # J/K
		amu2kg = 1.660539e-27 #amu/kg

		# constants with self.
		eo, ro = re.findall('van der Waals scaling parameters: eo= (.*?)eV, ro= (.*?)A', data)[0]
		self.ro = float(ro) # angstrom
		self.eo = float(eo) # eV
		self.piro2 = np.pi*(self.ro)**2 # Q = piro2 * Qst
		self.Tfac = kB/(self.eo*elC) # Tst = Tfac * T
		
		coll_mass = {'He': 4.0026, 'N2': 28.0062}
		
		self.direc = directory
		self.filename = re.findall('input file name = (.*?)\n',data)[0].strip()[:-4]
		
		#get atom mass
		try:
			self.mass = float(re.findall('mass of ion = (.*?)\n',data)[0].replace('D','E').replace('amu','').strip())
		except:
			self.mass = 0.00
		
		#get collision gas
		try:
			self.coll_gas = re.findall('Mobility Calculated under (.*?) gas',data)[0].strip()
		except:
			self.coll_gas = 'N2'
		
		#compute reduced mass:
		self.redmass = coll_mass[self.coll_gas]*self.mass / (coll_mass[self.coll_gas]+self.mass)
		self.gfac = np.sqrt( (self.eo*elC) / (0.5*self.redmass*amu2kg) ) # g = gfac * gst
		
		#get number of atoms
		try:
			self.Natoms = int(re.findall('number of atoms = (.*?)\n',data)[0].strip())
		except:
			self.Natoms = 0
		
		#get NTgrid
		try:
			self.NTgrid = int(re.findall('# of T_eff grid points = (.*?)\n',data)[0].strip())+1
		except:
			self.NTgrid = 1
		
		#get empirical correction data
		if 'Empirical correction turned off' in data:
			self.EmpCorr = False
			self.EmpCorr_A = 0.
			self.EmpCorr_B = 0.
		else:
			self.EmpCorr = True
			self.EmpCorr_A = float(re.findall('A =   (.*?)\n',data)[0].strip())
			self.EmpCorr_B = float(re.findall('B = (.*?) Td\n',data)[0].strip())
		
		#get integration parameters
		try:
			self.itn = int(re.findall('number of complete cycles (.*?)\n', data)[0].rsplit(' ')[-1].strip())
			self.inp = int(re.findall('number of velocity points (.*?)\n', data)[0].rsplit(' ')[-1].strip())
			self.imp = int(re.findall('number of random points (.*?)\n', data)[0].rsplit(' ')[-1].strip())
		except:
			#print('Problems reading itn, inp, imp')
			error_popup('critical','Value Error',f'Problems reading itn, inp, and imp - did {file} finish correctly?')
			return		 
		
		# momentum transfer data
		# gst q1st +/- q1st_err q2st +/- q2st_err q3st +/- q3st_err
		try:
			ifind = data.find('Final averaged values of Q*(l):')
			MT_str = ((data[ifind:ifind+112+self.inp*88]).replace('+/-','').split('\n')[3:])
			MT_arr = np.array([[float(s) for s in MT_str_line.split()] for MT_str_line in MT_str])
			if MT_arr.shape[0] != self.inp:
				print('not all MT data read in!')
			self.gst	  = MT_arr[:,0]	   # velocity grid
			self.qist	 = MT_arr[:,[1,3,5]] # Q^(l), l=1,2,3
			self.qist_err = MT_arr[:,[2,4,6]] # sig(Q^(l)), l=1,2,3
		except:
			error_popup('critical','Value Error',f'Problems reading momentum transfer data - did {file} finish correctly?')
			return

		# mobility summary
		# Teff [K]   E/N [Td]  K0 [cm**2/Vs]  CCS [A**2]  uncertainty
		try:
			mob_str = ((data[data.find('Mobility Summary')+25:]).replace('%','').split('\n')[4:-3])
			mob_arr = np.array([[float(s) for s in mob_str_line.split()] for mob_str_line in mob_str])
			self.Teff_r   = mob_arr[:,0]
			self.EN_r	 = mob_arr[:,1]
			self.K0_r	 = mob_arr[:,2]
			self.CCS_r	= mob_arr[:,3]
			self.uncert_r = mob_arr[:,4] # in percent!
			self.alpha_func = self.K0_r/self.K0_r[0] - 1.
			
			# check if there are more than just low-field data
			if len(self.Teff_r) > 1:
				self.highfield = True
			else:
				self.highfield = False
		except:
			#raise ValueError('Could not retrieve Mobility Summary. Job not finished!')
			error_popup('critical','Value Error',f'Problems reading Mobility summary table - did {file} finish correctly?')
			return
		
		# get calculation time
		try:
			time = float(re.findall('Job Completed in(.*?)\n',data)[0].strip().replace('s',''))
			self.time = str(round(time,2))
		except:
			self.time = 0.0
		
		# make a summary text to be printed in the GUI
		stext = []
		stext.append(self.filename.rsplit('.')[0])
		stext.append('')
		stext.append('itn=%i, inp=%i, imp=%i' %(self.itn, self.inp, self.imp))
		stext.append('Collision gas: %s' %self.coll_gas)
		stext.append('Empirical Correction: %s' %(self.EmpCorr))
		stext.append('')
		stext.append('T_bath set to %.0f K' %(self.Teff_r[0]))
		stext.append('CCS(T_bath) = %6.1f +/- %5.1f Ang^2' %(self.CCS_r[0], 1e-2*self.CCS_r[0]*self.uncert_r[0]))
		stext.append('K0( E/N=0 ) = %6.3f +/- %5.3f cm^2/Vs' %(self.K0_r[0], 1e-2*self.K0_r[0]*self.uncert_r[0]))
		stext.append('')
		if self.highfield:
			stext.append('CCS	  calculated up to %4.0f K' %(self.Teff_r[-1]))
			stext.append('Mobility calculated up to %4.1f Td' %(self.EN_r[-1]))
		
		self.summary_text = stext
	
	# function to calculate OM^(l,s) at Teff on the fly.
	def OM(self,l,s,Teff):
		'''Computes collision integral OM^(l,s) and its error at Teff (K) in Ang**2. l=1-3, s=1-4'''
		if Teff<self.Teff_r.min() or Teff>self.Teff_r.max():
			error_popup('critical','Value Error',f'Teff out of range! Must be between {self.Teff_r.min():.0f} and {self.Teff_r.max():.0f} K')
			return		 
		
		Tst = Teff*self.Tfac
		dgst = self.gst[1] - self.gst[0]
		w = 2./( factorial(s+1)*Tst**(s+2)) * self.gst**(2*s+3) * np.exp(-self.gst**2/Tst)
		q = self.qist[:,l-1]
		qerr = self.qist_err[:,l-1]
		# integrate
		I = np.sum(w*q*dgst) * self.piro2
		Ierr = np.sqrt(np.sum((w*qerr*dgst)**2)) * self.piro2
		
		return I, Ierr
	
	# function to calculate and plot the alpha function and coefficients
	# again, input check not necessary since order N is chosen from dropdown box
	def get_alpha_coeff(self,nord=6,plotting=False):
		'''Express K(E/N) = K(0)*[1 + a2*(E/N)^2 + a4*(E/N)^4 + ...]
		and get a2 in 1/Td^2 and a4 in 1/Td^4. nord is either 4 or 6'''
		K0_SI = self.K0_r * 1e-4 # in m^2/Vs
		
		# check order or approx (nord=4 or =6)
		if nord == 4:
			fit = lambda x, a2, a4: K0_SI[0]*(1. + a2*x**2 + a4*x**4)
		elif nord == 6:
			fit = lambda x, a2, a4, a6: K0_SI[0]*(1. + a2*x**2 + a4*x**4 + a6*x**6)
		else:
			error_popup('critical','Polynomial fitting Error','Problem with polynomial fitting; polynomial does not converge to the 4th or 6th order - check your mobility data in the .mout file and see if it makes sense.')
			return		 
					
		# fit to function
		popt, pcov = curve_fit(fit, self.EN_r, K0_SI)
		
		if plotting:
			# alpha coeff label for plot
			al_lab = ''
			for i in range(nord//2):
				a_exp = np.floor(np.log10(np.abs(popt[i]))) # round to next lower integer (also when negative!)
				a_pre = popt[i]/(10**a_exp) # prefactor of exponent
				al_lab += r'$\alpha_{%i}$ = $%+.4f\times10^{%i}$ Td$^{-%i}$'%(2*(i+1),a_pre,a_exp,2*(i+1))
				al_lab += '\n'
			
			# plot for visualization
			plt.figure()
			plt.plot(self.EN_r, self.K0_r, 'kx-', label='MobCal-MPI data')
			plt.plot(self.EN_r, 1e4*fit(self.EN_r, *popt), 'r--', label=al_lab)
			plt.legend()
			plt.xlabel(r'$E/N$ [Td]')
			plt.ylabel(r'$K_0$ [cm$^2$/Vs]')
			plt.tight_layout()
			plt.show(block=False)
		
		return popt
	
	## exporting functions
	def export_CCS(self):
		'''writes the CCS data to a file'''
		X = np.vstack([self.Teff_r, self.CCS_r, self.CCS_r*self.uncert_r*1e-2]).T
		
		try:
			fo = open(self.direc+self.filename+'_CCSdata.csv','w')
			fo.write('Teff [K],  CCS [A**2],  errCCS [A**2] \n')
			for Xi in X:
				fo.write(' %7.2f,   %7.2f,	  %5.2f' %(*Xi,) + '\n')
			fo.close()

		except(PermissionError):
			error_popup('warning','Permission error',f'Another instance of {self.filename}_CCSdata.csv is open. Please close it and retry exporting the summary.')
			return				
	
	def export_Ql(self):
		'''writes the momentum transfer data to a file'''
		g = self.gfac * self.gst
		Ql = self.qist * self.piro2
		sigQl = self.qist_err * self.piro2
		X = np.vstack([g, Ql[:,0], sigQl[:,0], Ql[:,1], sigQl[:,1], Ql[:,2], sigQl[:,2]]).T
		
		try:
			fo = open(self.direc+self.filename+'_Mom_Transfer_data.csv','w')
			fo.write('g [m/s],  Q(1) [A**2],  errQ(1) [A**2],  Q(2) [A**2],  errQ(2) [A**2],  Q(3) [A**2],  errQ(3) [A**2] \n')
			for Xi in X:
				fo.write('%6.1f,	 %6.2f,		%6.2f,		 %6.2f,		%6.2f,		 %6.2f,		%6.2f' %(*Xi,) + '\n')
			fo.close()

		except(PermissionError):
			error_popup('warning','Permission error',f'Another instance of {self.filename}_Mom_transfer_data.csv is open. Please close it and retry exporting the summary.')
			return   
					 
	def export_K(self):
		'''writes the mobility data to a file'''
		
		X = np.vstack([self.EN_r, self.K0_r, self.alpha_func]).T
		try:
			fo = open(self.direc+self.filename+'_Mobidata.csv','w')
			fo.write('E/N [Td],  K0 [cm2/Vs],  alpha func\n')
			for Xi in X:
				fo.write(' %6.2f,	 %6.4f,	  %+6.4f' %(*Xi,) + '\n')
			fo.close()

		except(PermissionError):
			error_popup('warning','Permission error',f'Another instance of {self.filename}_Mobidata.csv is open. Please close it and retry exporting the summary.')
			return		   
				 
	def export_summary(self):
		'''writes the summary table at the end of the .mout file to a new file'''
		X = np.vstack([self.Teff_r, self.EN_r, self.K0_r, self.CCS_r, self.uncert_r]).T
		
		try:
			fo = open(self.direc+self.filename+'_summary.csv','w')
			fo.write(' Teff [K],  E/N [Td], K0 [cm**2/Vs], CCS [A**2], uncertainty\n')
			for Xi in X:
				fo.write(' %7.2f,	%6.2f,	  %6.4f,	 %7.2f,	 %5.2f' %(*Xi,) + '%\n')
			fo.close()
		except(PermissionError):
			error_popup('warning','Permission error',f'Another instance of {self.filename}_summary.csv is open. Please close it and retry exporting the summary.')
			return				

	## purely plotting functions
	def plot_CCS(self):
		'''Plotting CCS over the available Teff range'''
		plt.figure(figsize=(4,3),dpi=250)
		plt.plot(self.Teff_r,self.CCS_r)
		plt.fill_between(self.Teff_r,self.CCS_r*(1.+self.uncert_r*1e-2),self.CCS_r*(1.-self.uncert_r*1e-2),color='C0',alpha=0.3)
		plt.xlabel(r'$T_{eff}$ [K]')
		plt.ylabel(r'CCS [$\AA^2$]')
		plt.tight_layout()
		plt.show(block=False)
	
	def plot_CCS_integrand(self):
		'''Plotting CCS integrand over the velocity range'''
		g = self.gfac * self.gst
		Tst = self.Teff_r[0]*self.Tfac # T_bath in T*
		l,s = 1,1
		w = 2./( factorial(s+1)*Tst**(s+2)) * self.gst**(2*s+3) * np.exp(-self.gst**2/Tst)
		Iy = self.qist[:,l-1]*w
		Iyerr = self.qist_err[:,l-1]*w
		plt.figure(figsize=(4,3),dpi=250)
		plt.plot(g,Iy,color='cornflowerblue',label=r'%.1f $\pm$ %.1f $\AA^2$' %(self.CCS_r[0], self.CCS_r[0]*self.uncert_r[0]*1e-2) )
		plt.fill_between(g,Iy+Iyerr,Iy-Iyerr,color='cornflowerblue',alpha=0.3)
		plt.xlabel(r'$g$ [m/s]')
		plt.ylabel(r'CCS Integrand [a.u.]')
		plt.legend()
		plt.yticks([])
		plt.tight_layout()
		plt.show(block=False)
	
	def plot_Qldat(self,lmax=3):
		'''Plots momentum transfer integrals Q^(l) up to lmax = 1,2 or 3.'''
		g = self.gfac * self.gst
		plt.figure(figsize=(4,3),dpi=250)
		for l in range(lmax): # l is actually l-1
			plt.plot(g,self.qist[:,l]*self.piro2,'C%i'%l,label=r'$l=%i$'%(l+1))
			plt.fill_between(g,(self.qist[:,l]+self.qist_err[:,l])*self.piro2,(self.qist[:,l]-self.qist_err[:,l])*self.piro2,color='C%i'%l,alpha=0.3)
		plt.ylabel(r'$Q^{(l)}$ [$\AA^2$]')
		plt.xlabel(r'$g$ [m/s]')
		plt.legend()
		plt.grid(which='both',lw=0.3)
		plt.tight_layout()
		plt.show(block=False)

class many_mout:
	def __init__(self,direc):
		self.direc = direc
		
		self.files = [x.rsplit('.',1)[0] for x in os.listdir(self.direc) if x.lower().endswith('.mout')]
		self.M_list = [mout_info(self.direc, file+'.mout') for file in self.files]
		self.Nfiles = len(self.files)
		self.maxlen_fname = np.max([len(file) for file in self.files])
		
		# make a summary text to be printed in the GUI
		stext = []
		stext.append('%i files read in' %self.Nfiles)
		
		# see if T_bath is the same
		Tbath_list = [int(M.Teff_r[0]) for M in self.M_list]

		if all(Tbath_list[0] == T for T in Tbath_list):
			stext.append('common Tbath: %i K' %Tbath_list[0])
			self.common_Tbath = True
		
		else:
			stext.append('WARNING! The %i .mout files contain jobs that have different bath gas temperatures that range from %.0f to %0.f K' %(self.Nfiles,np.min(Tbath_list),np.max(Tbath_list))+'. Please verify that this was intended before proceeding to analyze your results.')
			self.common_Tbath = False

			#warning message that bath gas temperatures are not consistent
			error_popup('warning','Temperature inconsistency',f'The {self.Nfiles} .mout files within the specified directory contain jobs that have different bath gas temperatures that range from {np.min(Tbath_list):.0f} to {np.max(Tbath_list):.0f} K. Please verify that this was intended before proceeding to analyze your results.')
			pass

		stext.append('')
		
		# check low/high field data
		if all([M.highfield==True for M in self.M_list]):
			stext.append('All files contain high field data')
			self.dat_type = 'highfield' # all mout files contain high field data
		elif all([M.highfield==False for M in self.M_list]):
			stext.append('All files contain low field data only')
			self.dat_type = 'lowfield' # all mout files contain only low field data
		else:
			stext.append('WARNING! The .mout directory contains files with high-field and exclusively low-field data!')
			stext.append('Consequently, only low field data can be exported!')
			self.dat_type = 'mixed' # some mout files have only low, some only high field data

			#print warning that mixed data is present
			error_popup('warning','Temperature inconsistency','The .mout files within the specified directory contain both low-field and high-field jobs! Only the low field data can be exported.')
			pass		 

		stext.append('')
		
		# all filenames with low field CCS and K0
		stext.append('filename' + ' '*(self.maxlen_fname-8) +' CCS [Ang**2]  K0 [cm2/Vs]')
		for M in self.M_list:
			stext.append("{0:<{1}}".format(M.filename, self.maxlen_fname) + '   %7.2f	   %6.4f' %(M.CCS_r[0],M.K0_r[0]))
		self.summary_text = stext
	
	def plot_CCS_list(self):
		'''plots all low field CCSs of M_list in a bar plot'''
		xpos = np.arange(self.Nfiles)
		CCSs = [M.CCS_r[0] for M in self.M_list]
		labels = [M.filename for M in self.M_list]
		
		plt.figure()
		h = plt.bar(xpos,CCSs)
		plt.subplots_adjust(bottom=0.3)
		xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in h]
		_ = plt.xticks(xticks_pos, labels, ha='right', rotation=45)
		plt.ylabel(r'CCS / $\AA^2$')
		plt.show(block=False)
	
	def plot_CCS_dist(self):
		'''plots the distribution (histogram) of all low field CCSs of M_list'''
		CCSs = [M.CCS_r[0] for M in self.M_list]
		minCCS = np.min(CCSs)
		maxCCS = np.max(CCSs)
		kde = gaussian_kde(CCSs)
		xCCS = np.linspace(minCCS*0.9,maxCCS*1.1,201)
		
		plt.figure()
		plt.hist(CCSs,bins=30,density=True,alpha=0.3,color='C0',label='histogram')
		plt.plot(xCCS,kde(xCCS),'C0-',label='Gaussian smooth')
		plt.xlabel(r'CCS / $\AA^2$')
		plt.ylabel('frequency')
		plt.legend(loc='upper right')
		plt.show(block=False)
	
	def export_CCS(self, itype=0):
		if self.common_Tbath:
			Tinfo = '# All CCS calculated at Tbath = %.1f K\n' % (self.M_list[0].Teff_r[0])
		else:
			Tinfo = '# WARNING: CCSs calculated at different Tbaths!\n'

		# Export low field data - mixed data can only be exported at lowfield. itype = 0 = lowfield, itype = 1 is for high-field
		if itype == 0:
			CCSs = [M.CCS_r[0] for M in self.M_list]
			errCCSs = [M.CCS_r[0] * M.uncert_r[0] * 1e-2 for M in self.M_list]

			fo = open(self.direc + 'Export_mout_CCS_lowfield.csv', 'w')
			fo.write('filename,T_bath[K],CCS [A**2],errCCS [A**2]\n')
			for i in range(self.Nfiles):
				fo.write("{0:<{1}}".format(self.files[i], self.maxlen_fname) + ',%7.2f,%7.2f,%5.2f' % (self.M_list[i].Teff_r[0], CCSs[i], errCCSs[i]) + '\n')

		# export high field data
		else:
			Teffs = [list(M.Teff_r) for M in self.M_list]

			# Because we had to convert the Teff to normal list for proper list comprehension, np is going to cause issues. We need to flatten the list in order for the interpolation
			flattened_Teffs = [temp for sublist in Teffs for temp in sublist]

			# check if all Teffs are compatible
			TF_set = [np.array_equal(Teff_i, Teffs[0]) for Teff_i in Teffs]

			if all(TF_set):
				Tgrid_info = '# All CCS calculated on the same Teff grid\n'
				Teff_grid = Teffs[0]
				Ngrid = len(Teff_grid)
				CCSs = np.array([interp1d(M.Teff_r, M.CCS_r)(Teff_grid) for M in self.M_list])
				
			else:
				Tgrid_info = '# WARNING: Teff grids not matching. Interpolation on a common range used! Teff values outside the common range will be assigned NaN\n'

				Ngrid = np.max([len(M.CCS_r) for M in self.M_list])
				common_Tmin = np.max(np.min(flattened_Teffs))
				common_Tmax = np.min(np.max(flattened_Teffs))
				Teff_grid = np.linspace(common_Tmin, common_Tmax, Ngrid)

				# Clip Teff_grid values to be within the interpolation range
				Teff_grid = np.clip(Teff_grid, np.min(flattened_Teffs), np.max(flattened_Teffs))

				# show a warning to the user
				error_popup('warning','Temperature inconsistency','The Teff grids do not match. Interpolation on a common range is being used when exporting high-field CCSs and/or mobilities! Teff values outside the common range will be skipped')
				pass

			# Initialize CCSs with empty strings before the loop
			CCSs = np.empty((self.Nfiles, Ngrid), dtype=object)
			CCSs.fill('')  # Fill the array with empty strings

			# Now we can proceed with the interpolation loop
			for i, M in enumerate(self.M_list):
				# Find indices of valid Teff values within the interpolation range of the current file
				valid_Teff_indices = [j for j, T in enumerate(Teff_grid) if np.min(M.Teff_r) <= T <= np.max(M.Teff_r)]

				# Interpolate CCS values for valid Teff values within the interpolation range
				for j in valid_Teff_indices:
					# Interpolate CCS for the j-th valid Teff value
					CCSs[i, j] = interp1d(M.Teff_r, M.CCS_r)(Teff_grid[j])

			#Now we can write to file
			try:
				fo = open(self.direc + 'Export_mout_CCS_highfield.csv', 'w')
				fo.write(Tinfo)
				fo.write(Tgrid_info)

				# files in columns
				fo.write(r'Teff(K) \ CCS(A^2) of' + ', %s' * self.Nfiles % (*self.files,) + '\n')
				for i in range(Ngrid):
					# Convert the list of strings to a comma-separated string with empty entries
					ccs_str = ', '.join(['%7.2f' %ccs if ccs != '' else ' ' for ccs in CCSs[:, i]]) 
					fo.write('%7.2f			  ' % Teff_grid[i] + ', ' + ccs_str + '\n')
				
				fo.close()

			except(PermissionError):
					error_popup('warning','Permission Error','Export_mout_CCS_highfield.csv is open in another window. Please close it and re-export your data.')
					return 
			
	def export_Mobility(self, itype=0):
		if self.common_Tbath:
			Tinfo = '# All mobility data calculated at Tbath = %.1f K\n' % (self.M_list[0].Teff_r[0])
		else:
			Tinfo = '# WARNING: mobility data calculated at different Tbaths!\n'

		# Export low field data - mixed data can only be exported at lowfield. itype = 0 = lowfield, itype = 1 is for high-field
		if itype == 0:
			K0s = [M.K0_r[0] for M in self.M_list]
			errK0s = [M.K0_r[0] * M.uncert_r[0] * 1e-2 for M in self.M_list]

			fo = open(self.direc + 'Export_mout_K0_lowfield.csv', 'w')
			fo.write('filename,T_bath[K],K0 [cm^2/Vs],errK0 [cm^2/Vs]\n')
			for i in range(self.Nfiles):
				fo.write("{0:<{1}}".format(self.files[i], self.maxlen_fname) + ',%6.4f,%6.4f,%6.4f' % (
					self.M_list[i].Teff_r[0], K0s[i], errK0s[i]) + '\n')

		# Export high field data
		else:
			ENs = [list(M.EN_r) for M in self.M_list]

			# Because we had to convert the EN to a normal list for proper list comprehension, np is going to cause issues. We need to flatten the list in order for the interpolation
			flattened_ENs = [temp for sublist in ENs for temp in sublist]
  
			# E/N grid will never the the same for two molecules because their mobilities will be different. So it is unnecessary to check if the E/N is computed on the same grid. Make sure the user knows this!
			error_popup('warning','Temperature Inconsistency','Interpolation on a common range is being used when exporting high-field mobilities even in Teff values are computing on the same grid! This is because conversion of Teff to E/N depends on instrinsic ion properties, and hence, will never be the same for two different analytes. Also note that E/N values outside the common range will be skipped.')
			pass 
			 
			Ngrid = np.max([len(M.CCS_r) for M in self.M_list])
			common_EN_min = np.max(np.min(flattened_ENs))
			common_EN_max = np.min(np.max(flattened_ENs))
			EN_grid = np.linspace(common_EN_min, common_EN_max, Ngrid)

			# Clip EN_grid values to be within the interpolation range
			EN_grid = np.clip(EN_grid, np.min(flattened_ENs), np.max(flattened_ENs))

			# Initialize K0s with the special value before the loop
			K0s = np.empty((self.Nfiles, Ngrid), dtype=object)
			K0s.fill('')  # Fill the array with empty strings

			# Now we can proceed with the interpolation loop
			for i, M in enumerate(self.M_list):
				# Find indices of valid EN values within the interpolation range of the current file
				valid_EN_indices = [j for j, EN in enumerate(EN_grid) if np.min(M.EN_r) <= EN <= np.max(M.EN_r)]

				# Interpolate K0 values for valid EN values within the interpolation range
				for j in valid_EN_indices:
					# Interpolate K0 for the j-th valid EN value
					K0s[i, j] = interp1d(M.EN_r, M.K0_r)(EN_grid[j])

			#Now we can write to file
			try:
				fo = open(self.direc + 'Export_mout_K0_highfield.csv', 'w')
				fo.write(Tinfo)

				# files in columns
				fo.write(r'E/N(Td) \ K0(cm^2/Vs) of' + ', %s' * self.Nfiles % (*self.files,) + '\n')
				for i in range(Ngrid):
					# Convert the list of strings to a comma-separated string with empty entries
					#k0_str = ', '.join([str(k0) if k0 != '' else ' ' for k0 in K0s[:, i]]) #old
					k0_str = ', '.join(['%6.4f'%k0 if k0 != '' else ' ' for k0 in K0s[:, i]]) #new
					fo.write('%6.2f' % EN_grid[i] + ' ' * (24 - 6) + ', ' + k0_str + '\n')
				
				fo.close()

			except (PermissionError):
				error_popup('warning','Permission Error','Export_mout_K0_highfield.csv is open in another window. Please close it and re-export your data.')
				return
			
