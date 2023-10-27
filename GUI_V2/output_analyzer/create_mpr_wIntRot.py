# -*- coding: utf-8 -*-
"""
Create mpr file including IntRot correction

Created on Sun Sep 18 12:48:55 2022

@author: alexa
"""
#%% load in libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'lines.linewidth': 2})
from scipy.linalg import block_diag
from scipy.optimize import minimize

#import sys
import os
#import re

# general constants
kB = 1.380658E-23 #J/K boltzmann constant
h_SI = 6.62607e-34 #Js, Planck constant
h_bar = h_SI/(2.*np.pi) # Js reduced Planck constant
c_SI = 2.99782e+8 #m/s, speed of light

# conversions
J2Eh = 2.2937123e+17 # 1 Joule in hartree
bohr2m = 5.29177211e-11 # m/bohr
bohr2A = bohr2m * 1.0e+10 # Angstrom/bohr
amu2kg = 1.660539e-27 # kg/amu
cm2in = 0.3937 # inch/cm
eigval2s2 = 1e20/amu2kg/J2Eh # Eh/Ang**2.amu = 2.62549964E+29 s**-2
piro2 = np.pi * 3.043**2 # in angstrom

# cutoffs for overlaps
cutoff_1 = 80.
cutoff_2 = 40.
cutoff_12 = 60.

# colors of atoms
at_cols = {'C': 'black', 'H': 'dimgrey', 'O': 'orangered', 'N': 'royalblue',
           'S': 'gold', 'F': 'cyan', 'Cl': 'lime', 'Br': 'firebrick'}

# perturbation symbols e_ijk, Levi-Civita symbols in 3D
eps_LC = np.zeros((3,3,3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            if i == j or i == k or j == k:
                eps_LC[i,j,k] = 0
            elif [i,j,k] == [0,1,2] or [i,j,k] == [2,0,1] or [i,j,k] == [1,2,0]:
                eps_LC[i,j,k] = 1
            elif [i,j,k] == [0,2,1] or [i,j,k] == [1,0,2] or [i,j,k] == [2,1,0]:
                eps_LC[i,j,k] = -1

#%% misc functions
def print_mat(A,fmt='%+.3f'):
    '''Prints the matrix A in format fmt'''
    ffmt = fmt+' '
    #print('')
    for i in range(A.shape[0]):
        print(ffmt*A.shape[1] %(*A[i],))

def get_CoM(xyz,mass):
    ''' return the CoM of xyz (3xNat array)'''
    return np.sum(mass*xyz,axis=1)/np.sum(mass)

def symm_mat(M):
    '''symmetrizes a quadratic matrix M by taking 1/2(M[i,j]+M[j,i])'''
    return 0.5*(M + M.T)

def mat_sqrt(A,thresh=1e-5):
    '''returns matrix B such that B*B=A. Thus, B=A^1/2. A is assumed to be hermetian.'''
    N = A.shape[0]
    LAMB, V = np.linalg.eigh(A)
    LAMBr = np.real(LAMB)
    V = np.real(V)
    LAMBsq_arr = np.zeros(N)
    for i in range(N):
        if LAMBr[i] > thresh:
            LAMBsq_arr[i] = np.sqrt(LAMBr[i]) # keep order!
        else:
            LAMBsq_arr[i] = 0.
    LAMBsq = block_diag(*LAMBsq_arr)
    sqrtA = np.dot(V, np.dot(LAMBsq,V.T))
    return sqrtA

def mat_inv(A,thresh=1e-6):
    '''returns generalized inverse of singular matrix. A is assumed to be hermetian.'''
    N = A.shape[0]
    LAMB, V = np.linalg.eigh(A)
    LAMBr = np.real(LAMB)
    V = np.real(V)
    LAMBi_arr = np.zeros(N)
    for i in range(N):
        if np.abs(LAMBr[i]) > thresh:
            LAMBi_arr[i] = 1./LAMBr[i]
        else:
            LAMBi_arr[i] = 0.
    LAMBi = block_diag(*LAMBi_arr)
    Ai = np.dot(V, np.dot(LAMBi,V.T))
    return Ai

#%% molecule class
class molecule:
    def __init__(self,directory,basen,EHigh=True,Mobil=True,IntRot=False):
        '''Using basen.inp, basen.out, basen.hess, basen_IntHessPrint_internal.hess, and
        basen_IntHessPrint.opt, all relevant data are read it'''
        self.directory = directory
        self.basen = basen
        print('\n')
        print('== Reading in data from %s%s files ==' %(directory,basen))
        ### get charge, multiplicity, and energy from _OptFreq.out file
        file = basen + '_OptFreq.out'
        f = open(directory+file)
        lines = f.readlines()
        f.close()
        
        self.dipole = None # only read in for solvents
        self.polariz = None
        self.dipole_ax = None
        for line in lines:
            if line.startswith(' Total Charge'):
                self.charge = float(line.split()[-1])
            if line.startswith(' Multiplicity'):
                self.multi = float(line.split()[-1])
            if line.startswith('FINAL SINGLE POINT ENERGY'):
                self.Eelec = float(line.split()[-1])
            if line.startswith('Rotational constants in cm-1'):    #get rotational constants
                self.RotA = float(line.split()[4]) * 100.*c_SI*1E-9  # in GHz
                self.RotB = float(line.split()[5]) * 100.*c_SI*1E-9
                self.RotC = float(line.split()[6]) * 100.*c_SI*1E-9
            if line.startswith('Point Group:'): #get rotational symmetrie number
                self.sigma_OR = int(line.split()[-1]) # integer
            if line.startswith('Magnitude (Debye)'):
                self.dipole = float(line.split()[-1]) # in Debye
            if line.startswith('Isotropic polarizability :'):
                self.polariz = float(line.split()[-1])*1.4818e-31 # in m^3
            if line.startswith('x,y,z [Debye]:'):
                mu_abc = np.array([float(s) for s in line.split()[-3:]])
                self.dipole_ax = ['A','B','C'][np.argmax(np.abs(mu_abc))]
        
        ### get Hessian, full normal modes, geometry and masses
        file = basen + '_OptFreq.hess'
        f = open(directory+file)
        lines = f.readlines()
        f.close()
        
        # find blocks with different information
        for i in range(len(lines)):
            if lines[i].startswith('$hessian'):
                i_hess = i
            if lines[i].startswith('$atoms'):
                i_atom = i
            if lines[i].startswith('$normal_modes'):
                i_NM = i
            if lines[i].startswith('$vibrational_frequencies'):
                i_vibs = i
        
        self.Nat = int(lines[i_atom+1]) # get number of atoms
        self.DoV = 3*self.Nat-6 # number of internal motions, linear molecules not supported
        
        ## read atoms and masses
        atoms = np.zeros((self.Nat,3))
        self.atom_labs = []
        self.mass = np.zeros(self.Nat)
        for i in range(self.Nat):
            DAT = lines[i_atom+2+i].split()
            self.atom_labs.append(DAT[0])
            self.mass[i] = float(DAT[1]) # data is in amu
            atoms[i,0] = float(DAT[2]) # data is in Bohr
            atoms[i,1] = float(DAT[3])
            atoms[i,2] = float(DAT[4])
        atoms = atoms.T * bohr2A # convert to Angstrom
        
        print('\nOutput and Hessian file successfully read in!')
        # CoM and pverall inertia tensor
        CoM = get_CoM(atoms,self.mass)
        print('\nCoM = %+.6f %+.6f %+.6f' %(*CoM,))
        self.atoms = np.array([atoms[:,i]-CoM for i in range(self.Nat)]).T # remove CoM
        self.m_tot = np.sum(self.mass)
        self.II = np.zeros((3,3))
        for i in range(3):
            self.II[i,i] = np.sum([ self.mass[a]*( np.sum(self.atoms[:,a]**2) - self.atoms[i,a]**2) for a in range(self.Nat)])
        for i,j in [[1,0],[2,0],[2,1]]:
            self.II[i,j] = -np.sum([ self.mass[a]*self.atoms[i,a]*self.atoms[j,a] for a in range(self.Nat) ])
            self.II[j,i] = self.II[i,j]
        print('Inertia tensor (amu*A**2)')
        print_mat(self.II,fmt='%+8.3f')
        
        # mass matrices
        self.sqM = np.zeros((3*self.Nat,3*self.Nat)) # diagonal matrix with sqrt(m) as values
        self.sqMi = np.zeros((3*self.Nat,3*self.Nat)) # diagonal matrix with 1/sqrt(m) as values
        for i in range(3*self.Nat):
            sqm = np.sqrt(self.mass[i//3])
            self.sqM[i,i] = sqm
            self.sqMi[i,i] = 1./sqm
        
        ## read Hessian
        Hess = np.zeros((3*self.Nat,3*self.Nat))
        Nblox = (3*self.Nat)//5 + int((3*self.Nat)%5 != 0)
        for N in range(Nblox):
            for i in range(3*self.Nat):
                ind = lines[i_hess+2+N*(3*self.Nat+1)].split()
                data = lines[i+i_hess+3+N*(3*self.Nat+1)].split()
                for jj in range(len(ind)):
                    j = int(ind[jj])
                    Hess[i,j] = float(data[jj+1])
        
        self.Hess = symm_mat(Hess) # symmetrize Hessian
        
        ## read L matrix (dimensionless)
        self.Lfull = np.zeros((3*self.Nat,3*self.Nat))
        Nblox = (3*self.Nat)//5 + int((3*self.Nat)%5 != 0)
        for N in range(Nblox):
            for i in range(3*self.Nat):
                ind = lines[i_NM+2+N*(3*self.Nat+1)].split()
                data = lines[i+i_NM+3+N*(3*self.Nat+1)].split()
                for jj in range(len(ind)):
                    j = int(ind[jj])
                    self.Lfull[i,j] = float(data[jj+1])
        self.L = self.Lfull[:,6:]
        
        ## read original vibs from ORCA
        self.vibs_orig = np.zeros(self.DoV)
        for i in range(self.DoV):
            self.vibs_orig[i] = float(lines[i_vibs+8+i].split()[-1])
        self.w_pr = np.copy(self.vibs_orig) # in case we don't do any IntRot projection
        
        ##### optional data #####
        
        # higher level SPE calc
        if EHigh:
            print('\n= Updating electronic energy =')
            self.Eelec = self.update_Ehigh()
        
        # mobility data
        if Mobil:
            print('\n= Reading in mobility data =')
            self.QDAT = self.read_Qast()
            self.plot_Qldat()
            plt.close()
        
        # internal rotation correction
        self.N_IR = 0 # in case we are not looking for internal rotations
        if IntRot:
            print('\n= Reading in fragment information =')
            # read in internal hessian and Bmat
            self.RedDic, self.RedList, self.Bmat, self.Lq, self.Hq_read = self.read_IntHess()
            self.NRed = len(self.RedList)
            self.Nbond = len(self.RedDic["bonds"])
            self.Nangle = len(self.RedDic["angles"])
            self.Ndihed = len(self.RedDic["dihedrals"])
            self.Ntors = len(self.RedDic["torsions"])
            self.NCart = len(self.RedDic["cart"])
            # read fragement data and built fragment connectivity list
            self.Nfrag, self.frag_conn, self.rot_bonds, self.frag_ats = self.read_frag()
            frag_conx = []
            for i in range(self.Nfrag):
                frag_conx.append([i+1])
            for i in range(self.Nfrag-1):
                a,b = np.sort(self.frag_conn[i])
                frag_conx[a-1].append(b)
                frag_conx[b-1].append(a)
            # perform projection
            self.N_IR, self.w_pr, self.w_IR, self.Vec_pr, self.D_IR, self.sig_ir = self.IntRot_proj(self.rot_bonds)
            # get Iir and Bir
            self.Iir = np.zeros(self.N_IR)
            self.Bir = np.zeros(self.N_IR)
            # for each bond in rot_bonds,
            # find all fragments and respective atoms rotating against each other
            for i in range(self.Nfrag-1):
                F1_i, F2_i = self.frag_conn[i]
                
                # find all atoms rotating against each other
                self.output = [F1_i]
                self.get_connection(F1_i,F2_i,frag_conx)
                F1_list = np.copy(self.output)
                F1_ats = np.hstack([ self.frag_ats[j-1] for j in F1_list ])
                
                self.output = [F2_i]
                self.get_connection(F2_i,F1_i,frag_conx)
                F2_list = np.copy(self.output)
                F2_ats = np.hstack([ self.frag_ats[j-1] for j in F2_list ])
                
                F1_ats = np.hstack(F1_ats)
                F2_ats = np.hstack(F2_ats)
                print('\n= Fragments rotation: ', F1_list[::-1],'--',F2_list,' =')
                #self.Iir[i], self.sig_ir[i], self.Bir[i] = self.IntRot_eval(self.rot_bonds[i],F1_ats,self.sig_ir[i])
                self.Iir[i], self.sig_ir[i], self.Bir[i] = self.IntRot_eval(self.rot_bonds[i],F2_ats,self.sig_ir[i])
        
        #############----- reading done -----#############
        self.ZPE = 0.5*h_SI*c_SI*100.*np.sum(self.w_pr)*J2Eh
        
        #### __init__ function done ####
    
    def get_connection(self,F1,F2,frag_conx):
        '''get all fragments connecting to fragment F1 but not towards fragment F2 recursively'''
        for i in frag_conx[F1-1]:
            if i == F1 or i == F2:
                continue
            else:
                if i in self.output:
                    continue
                else:
                    self.output.append(i)
                    self.get_connection(i,F2,frag_conx)
            return
    
    def plot_NMD_3D(self,L,extra=np.zeros((3,2)),col='b',cov_thresh=1.55):
        '''plots normal mode dispacement L in Cartesian coordinates'''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        # atoms and bonds
        for a in range(self.Nat):
            ax.scatter(self.atoms[0,a],self.atoms[1,a],self.atoms[2,a],color=at_cols[self.atom_labs[a]],
                       edgecolor='k',alpha=1.0,s=self.mass[a]*5)
        ax.set_box_aspect((np.ptp(self.atoms[0]), np.ptp(self.atoms[1]), np.ptp(self.atoms[2])))
        for bond in self.RedDic["bonds"]:
            a1,a2 = bond
            p1 = self.atoms[:,a1]
            p2 = self.atoms[:,a2]
            pm = p1+0.5*(p2-p1) # half way between atoms a1 and a2
            if np.linalg.norm(p2-p1) < cov_thresh: # covalent bonds
                ax.plot([p1[0],pm[0]],[p1[1],pm[1]],[p1[2],pm[2]],at_cols[self.atom_labs[a1]],alpha=0.8)
                ax.plot([pm[0],p2[0]],[pm[1],p2[1]],[pm[2],p2[2]],at_cols[self.atom_labs[a2]],alpha=0.8)
            else: # non-covalent bonds
                ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],'k--', alpha=0.5)
        # NMD
        for a in range(self.Nat):
            vec = L[3*a:3*(a+1)]
            xyz_s = self.atoms[:,a]
            xyz_e = xyz_s + vec
            ax.plot([xyz_s[0],xyz_e[0]],[xyz_s[1],xyz_e[1]],[xyz_s[2],xyz_e[2]],'r')
        ax.plot([xyz_s[0],xyz_e[0]],[xyz_s[1],xyz_e[1]],[xyz_s[2],xyz_e[2]],'r',label='displacement')
        # extra data
        ax.plot(extra[0],extra[1],extra[2],color=col,alpha=0.8,label='rotating axis/bond')
        # CoM
        ax.scatter([0.],[0.],[0.],marker='x',color='k',s=20,label='CoM')
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.legend()
        plt.show()
    
    def update_Ehigh(self):
        '''updates the electronic energy from a higher level calculation'''
        f = open(self.directory + self.basen + '_Ehigh.out', 'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if line.startswith('FINAL SINGLE POINT ENERGY'):
                Ehigh = float(line.split()[-1])
        
        print('E_high = %.10f Eh' %Ehigh)
        return Ehigh

    def read_Qast(self):
        '''Function to read in mobility data from MobCal-MPI 2.0 *.mout file. fname is without .mout extension'''
        f = open(self.directory+self.basen+'.mout')
        LINES = f.readlines()
        f.close()
        
        # get inp and section where Q* data is printed
        i_find = 0
        for i,line in enumerate(LINES):
            if line.startswith(' number of velocity points'):
                inp = int(line.split()[-1])
            if 'Final averaged values of Q*(l):' in line:
                i_find = i
    
        # read in momentum transfer cross sections on g* grid including errors
        gst = []
        q1st = []
        q2st = []
        q3st = []
        q1st_err = []
        q2st_err = []
        q3st_err = []
        for i in range(inp):
            lline = LINES[i_find+3+i].split()
            gst.append(float(lline[0]))
            q1st.append(float(lline[1]))
            q1st_err.append(float(lline[3]))
            q2st.append(float(lline[4]))
            q2st_err.append(float(lline[6]))
            q3st.append(float(lline[7]))
            q3st_err.append(float(lline[9]))
        
        QLDATA = np.array([gst,q1st,q2st,q3st,q1st_err,q2st_err,q3st_err]).T
        
        return QLDATA
    
    def plot_Qldat(self):
        # plotting
        gst,q1st,q2st,q3st,q1st_err,q2st_err,q3st_err = self.QDAT.T
        plt.figure(figsize=(10,8))
        plt.loglog(gst,q1st*piro2,'C0',label=r'$l=1$')
        plt.loglog(gst,q2st*piro2,'C1',label=r'$l=2$')
        plt.loglog(gst,q3st*piro2,'C2',label=r'$l=3$')
        plt.fill_between(gst,(q1st+q1st_err)*piro2,(q1st-q1st_err)*piro2,color='C0',alpha=0.3)
        plt.fill_between(gst,(q2st+q2st_err)*piro2,(q2st-q2st_err)*piro2,color='C1',alpha=0.3)
        plt.fill_between(gst,(q3st+q3st_err)*piro2,(q3st-q3st_err)*piro2,color='C2',alpha=0.3)
        plt.ylabel(r'$Q^{(l)}$ / $\AA^2$')
        plt.xlabel(r'$g^\ast$')
        plt.legend()
        plt.grid(which='both',lw=0.3)
        plt.tight_layout()
        plt.savefig(self.directory+self.basen+'_Ql.png', dpi=400,
                    format='png', bbox_inches='tight')
        print('Momentum Transfer Cross Section plot saved to %s' %self.directory)
        plt.show()
        
    def read_IntHess(self):
        '''Read in internal hessian and B matrix from basen_IntHessPrint.opt and
        basen_IntHessPrint_internal.hess files'''
        ### read internal coordinates and B-matrix
        file = self.basen + '_IntHessPrint.opt'
        f = open(self.directory+file)
        lines = f.readlines()
        f.close()
        
        for i in range(len(lines)):
            if lines[i].startswith('$redundant_internals'):
                i_RedList = i
            if lines[i].startswith('$bmatrix'):
                i_BMat = i
                break
        
        ## get dimensions
        NRed = int(lines[i_BMat+1].split()[0]) # number of redundant internals
        Nbond, Nangle, Ndihed, Ntors, NCart = [int(s) for s in lines[i_RedList+1].split()]
        
        ## read the list of internals
        RedDic = {"bonds": [], "angles": [], "dihedrals": [], "torsions": [], "cart": []}
        RedList = []
        for i in range(Nbond):
            Bond = [ int(s) for s in lines[i_RedList+2+i].split()[:-1] ]
            RedDic["bonds"].append(Bond)
            RedList.append(Bond)
        for i in range(Nangle):
            Angle = [ int(s) for s in lines[i_RedList+2+Nbond+i].split()[:-2] ]
            RedDic["angles"].append(Angle)
            RedList.append(Angle)
        for i in range(Ndihed):
            Dihed = [ int(s) for s in lines[i_RedList+2+Nbond+Nangle+i].split()[:-1] ]
            RedDic["dihedrals"].append(Dihed)
            RedList.append(Dihed)
        for i in range(Ntors):
            Torsion = [ int(s) for s in lines[i_RedList+2+Nbond+Nangle+Ndihed+i].split()[:-1] ]
            RedDic["torsions"].append(Torsion)
            RedList.append(Torsion)
        for i in range(NCart):
            line = lines[i_RedList+2+Nbond+Nangle+Ndihed+Ntors+i].split()
            Cart = [ int(line[0]), line[1] ] # [1, 'x']
            RedDic["cart"].append(Cart)
            RedList.append(Cart)
        
        ## read B-matrix
        Bmat = np.zeros((NRed,3*self.Nat))
        Nblox = (3*self.Nat)//6+int((3*self.Nat)%6 != 0)
        for N in range(Nblox):
            ind = [int(s) for s in lines[i_BMat+2+N*(NRed+1)].split()]
            for i in range(NRed):
                data = lines[i+i_BMat+3+N*(NRed+1)].split()
                for jj in range(len(ind)):
                    j = ind[jj]
                    Bmat[i,j] = data[jj+1]
        
        # mw-normal modes in internal coordinates
        Lqfull = np.dot(Bmat, np.dot(self.sqM,self.Lfull))
        #Lqfull = np.dot(Bmat, Lfull)
        Lq = Lqfull[:,6:] # first 6 modes are zero
        
        # normalize Lq matrix
        for i in range(self.DoV):
            Lq[:,i] /= np.linalg.norm(Lq[:,i])
        
        ## read in Hessian in internal coordinates, Hq
        file = self.basen + '_IntHessPrint_internal.hess'
        f = open(self.directory+file)
        lines = f.readlines()
        f.close()
        
        for i in range(len(lines)):
            if lines[i].startswith('$hessian'):
                i_hess_int = i
        
        Hq_read = np.zeros((NRed,NRed))
        Nblox = int(NRed/5)+int((NRed)%5 != 0)
        for N in range(Nblox):
            for i in range(NRed):
                ind = lines[i_hess_int+2+N*(NRed+1)].split()
                data = lines[i+i_hess_int+3+N*(NRed+1)].split()
                for jj in range(len(ind)):
                    j = int(ind[jj])
                    Hq_read[i,j] = float(data[jj+1])
        
        return RedDic, RedList, Bmat, Lq, Hq_read
    
    def read_frag(self):
        '''Reads fragment information from basen_IntHessPrint.inp file'''
        file = self.basen + '_IntHessPrint.inp'
        f = open(self.directory+file)
        lines = f.readlines()
        f.close()
        
        frag_flag = False
        for i,line in enumerate(lines):
            if 'ConnectFragments' in line:
                i_frag = i
                frag_flag = True
            if '* xyz' in line:
                i_xyz = i
        
        Nfrag = 1
        if frag_flag: # only if molecule actually has fragments!
            i=1
            line = lines[i_frag+i]
            frag_conn = []
            rot_ats = []
            while not ('end' in line):
                dat = line.split()
                frag_conn.append( [int(dat[1]),int(dat[2])] )
                rot_ats.append( [int(dat[4]),int(dat[5])] )
                Nfrag = np.max([ Nfrag, np.max([int(dat[1]),int(dat[2])]) ])
                i+=1
                line = lines[i_frag+i]
            print('%i fragments found!' %Nfrag)
            
            frag_ats = []
            for i in range(Nfrag):
                frag_ats.append([])
            for i in range(self.Nat):
                dat = lines[i_xyz+i+1]#.split()[0]
                frag_i = int( dat[dat.find('(')+1:dat.rfind(')')] )
                frag_ats[frag_i-1].append(i) # append current atom number to fragment list
        
        return Nfrag, frag_conn, rot_ats, frag_ats

    def IntRot_proj(self,rot_search_list):
        '''project out internal rotations and return updated vibs'''
        ### general checking of internal rotations
        ## method 1 of [Ayala1998]
        print('\n== Internal Rotation Projection==')
        print('\nMethod 1:')
        print('cutoff: %.0f' %cutoff_1 + '%')
        for i in range(self.Nbond): # loop over all bonds in RedDir
            bond_ats = self.RedDic["bonds"][i]
            for r in range(self.DoV): # loop over all modes
                T2 = 0.
                for j in range(self.Ndihed): # loop over all dihedrals
                    rot_ats = self.RedDic["dihedrals"][j][1:3]
                    if bond_ats == rot_ats or bond_ats == rot_ats[::-1]:
                        T2 += self.Lq[self.Nbond+self.Nangle+j,r]**2
                T = np.sqrt(T2)*100.
                if T > cutoff_1:
                    print('mode %2i is composed of %.1f' %(r,T)
                          + "% " + "D(A,%2i,%2i,B) dihedral changes" %(*bond_ats,))
        
        ## method 2 of [Ayala1998]
        # G, Gi and projector P=GGi
        u = self.sqMi**2
        G = np.dot(self.Bmat, np.dot(u,self.Bmat.T))
        Gi = mat_inv(G)
        P = np.dot(G,Gi)
        
        ## remove redundancy in coordinates
        #G_val, G_vec = np.linalg.eigh(G)
        #nonzero_G_val = [i for i in range(NRed) if G_val[i]>1e-7]
        #Kmat = G_vec[:,nonzero_G_val]
        #B_nr = np.dot(Kmat.T,Bmat) # non-redundant B-matrix of shape DoVx(3Nat)
        #Bi_nr = np.dot(B_nr.T,mat_inv(np.dot(B_nr,B_nr.T)))
        ## np.dot(B_nr,Bi_nr) is unity matrix
        
        Bi = np.dot(np.dot(u,self.Bmat.T),Gi) # this is Bi.T in [Ayala1998]
        Hq = np.dot(Bi.T,np.dot(self.Hess,Bi)) # diff to Hq_read is small-ish (good!)
        #Hq = np.copy(Hq_read)
        
        if self.NRed == self.DoV:
            P = np.identity(self.NRed)
        
        # constrain matrix
        if rot_search_list == 'all':    # all dihedrals
            C_dihed = np.zeros(self.Ndihed)
        else:                           # only dihedrals in rot_search_list
            C_dihed = np.ones(self.Ndihed)
            print('\nUnconstrained dihedrals:')
            for i in range(self.Ndihed):
                bond = self.RedDic["dihedrals"][i][1:3]
                if bond in rot_search_list or bond[::-1] in rot_search_list:
                    C_dihed[i] = 0.
                    print(i, self.RedDic["dihedrals"][i])

        C_arr = np.hstack(([1.]*self.Nbond,[1.]*self.Nangle,C_dihed,[1.]*self.Ntors,[1.]*self.NCart))
        C = block_diag(*C_arr)
        
        # new projector, P'
        CPC = np.dot(C, np.dot(P,C))
        CPCi = mat_inv(CPC)
        PC = np.dot(P,C)
        CP = np.dot(C,P)
        P_pr = P - np.dot(PC, np.dot(CPCi,CP))
        
        # new G, G'
        CGC = np.dot(C, np.dot(G,C))
        CGCi = mat_inv(CGC)
        GC = np.dot(G,C)
        CG = np.dot(C,G)
        G_pr = G - np.dot(GC, np.dot(CGCi,CG))
        
        # transform Hessian
        Hq_pr = np.dot(P_pr,np.dot(Hq,P_pr))
        #GH = np.dot(G_pr,Hq_pr)
        sqG_pr = mat_sqrt( G_pr )
        
        GHG = np.dot(sqG_pr,np.dot(Hq_pr,sqG_pr)) # mass weighted projected Hessian
        
        lamb_o, VEC = np.linalg.eigh(GHG) # VECs are all normal modes when only dihedral changes are allowed
        sort_IR = np.argsort(np.abs(lamb_o))
        lamb_r = lamb_o[sort_IR] / bohr2A**2
        VEC_IR_full = VEC[:,sort_IR]
        
        lamb = np.round(lamb_r,8)*eigval2s2
        w_IR_full = np.sign(lamb)*np.sqrt(np.abs(lamb))/(2.*np.pi*c_SI*100.) # in cm**-1
        
        Rot_vec_list = [i for i in range(self.NRed) if np.abs(w_IR_full[i])>2.] # which are actually >0
        N_IR = len(Rot_vec_list) # number of internal rotations
        w_IR = w_IR_full[Rot_vec_list] # vibrational freqs of internal rotations
        VEC_IR = VEC_IR_full[:,Rot_vec_list] # mw-normal coordinates of int-rots
        gamma_IR = 2.*np.pi*c_SI*100.*w_IR/(h_bar) *amu2kg*1e-20 # in 1/amu.A**2
        
        # print info of Method 2
        print('\nMethod 2:')
        print("rank of G^1/2 Hq' G^1/2  is  %2i  (equal to rotatable bonds)" %np.linalg.matrix_rank(GHG))
        print('\ncutoff: %.0f' %cutoff_2 + '%')
        print('Int Rot    Freq/cm**-1       Vibration    Freq/cm**-1     Overlap/%')
        print('---------------------------------------------------------------------')
        for ii in range(N_IR):
            for j in range(self.DoV):
                overlap = np.dot(VEC_IR[:,ii],self.Lq[:,j])*100.
                if np.abs(overlap) > cutoff_2:
                    print(' %2i        %7.2f              %2i          %7.2f          %5.1f' %(ii,w_IR[ii],j,self.vibs_orig[j],np.abs(overlap)))
        
        ## method 1+2 of [Ayala1998]
        print('\nMethod 1+2:')
        print('cutoff: %.0f' %cutoff_12 + '%')
        bond_ats_save = [0]*N_IR
        T_save = [0]*N_IR
        for i in range(self.Nbond): # loop over all bonds in RedDir
            bond_ats = self.RedDic["bonds"][i]
            for r in range(N_IR): # loop over all internal rotations
                T2 = 0.
                for j in range(self.Ndihed): # loop over all dihedrals
                    rot_ats = self.RedDic["dihedrals"][j][1:3]
                    if bond_ats == rot_ats or bond_ats == rot_ats[::-1]:
                        T2 += VEC_IR[self.Nbond+self.Nangle+j,r]**2
                T = np.sqrt(T2)*100.
                if T > cutoff_12:
                    print('IntRot mode %2i is composed of %.1f' %(r,T)
                          + "% " + "D(A,%2i,%2i,B) dihedral changes" %(*bond_ats,))
                if T > T_save[r]:
                    bond_ats_save[r] = bond_ats
                    T_save[r] = T
        
        print('\nRot vec list:')
        print(Rot_vec_list)
        
        
        ###############################################################
        ## projection of overall trans, overall rot and internal rot ##
        ## using projection matrix defined by [Szalay2014]           ##
        ###############################################################
        
        ## basis vectors T,R,P,Q of full space
        Dvar = np.zeros((3*self.Nat,3*self.Nat))
        
        ## for all IntRots, we need the da/dq matrix, where we choose q to be the VEC_IR coordinates
        ## da/dq must be w/o mass weighting since that is done within the eq. for the D vectors!
        ## VEC_IR is mass weighted so we need to remove that with np.dot(sqG_pr,VEC_IR)
        ## then renormalize the vectors
        ## backtransformation via Bi is non-linear but might be sufficient here
        Lq_IR = np.dot(sqG_pr,VEC_IR)
        for i in range(N_IR):
            Lq_IR[:,i] /= np.linalg.norm(Lq_IR[:,i])
        
        Bi_var = np.dot(Bi, Lq_IR )
        #Bi_var = np.dot(sqM,np.dot(Bi, Lq_IR ) )
        
        ## normalize Bi_var
        #for i in range(N_IR):
        #    Bi_var[:,i] /= np.linalg.norm(Bi_var[:,i])
        #    Bi_var[:,i] /= np.sqrt(np.abs(gamma_IR[i]))
        
        #Bi_var = np.dot(sqMi,Bi_var)
        
        #### manual overwriting
        #Bi_var = np.copy(L[:,0:1])
        #N_IR = 1
        #w_IR = np.array(vibs[0:1])
        ####
        
        ## generalized 3+N_IR dimensional inertia tensor
        I_gen = np.zeros((3+N_IR,3+N_IR))
        
        # inertia tensor entries for ax-ax (eq. 11a)
        I_gen[:3,:3] = np.copy(self.II)
        
        # inertia tensor entries for ax-int (eq. 11b) -> sometimes they vanish
        for i in range(3):
            for j in range(N_IR):
                for i1 in range(3):
                    for i2 in range(3):
                        I_gen[i,3+j] += eps_LC[i,i1,i2]*np.sum([ self.mass[at]*self.atoms[i1,at]*Bi_var[i2+3*at,j] for at in range(self.Nat) ])
                I_gen[3+j,i] = I_gen[i,3+j]
        
        # inertia tensor entries for int-int (eq. 11c)
        for i in range(N_IR):
            for j in range(N_IR):
                for al in range(3):
                    I_gen[3+i,3+j] += np.sum([ self.mass[at]*Bi_var[al+3*at,i]*Bi_var[al+3*at,j] for at in range(self.Nat) ])
        
        # form inverse-squareroots of generalized inertia tensor
        #I_gen_is = mat_sqrt(np.linalg.inv(I_gen))
        I_gen_is = np.linalg.inv(mat_sqrt(I_gen,thresh=1e-8))
        
        ## overall Trans (eq. 12a)
        Dtrans = np.zeros((3*self.Nat,3))
        M_sqrt = np.sqrt(self.m_tot)
        Dtrans[:,0] = np.hstack([[np.sqrt(self.mass[a]),0.,0.] for a in range(self.Nat)])/M_sqrt
        Dtrans[:,1] = np.hstack([[0.,np.sqrt(self.mass[a]),0.] for a in range(self.Nat)])/M_sqrt
        Dtrans[:,2] = np.hstack([[0.,0.,np.sqrt(self.mass[a])] for a in range(self.Nat)])/M_sqrt
        
        ## overall Rotation (eq. 12b)
        Drot = np.zeros((3*self.Nat,3))
        for al in range(3):
            for at in range(self.Nat):
                # ax-ax terms
                Drot[3*at+0,al] = (I_gen_is[al,1]*self.atoms[2,at]-I_gen_is[al,2]*self.atoms[1,at])*np.sqrt(self.mass[at])
                Drot[3*at+1,al] = (I_gen_is[al,2]*self.atoms[0,at]-I_gen_is[al,0]*self.atoms[2,at])*np.sqrt(self.mass[at])
                Drot[3*at+2,al] = (I_gen_is[al,0]*self.atoms[1,at]-I_gen_is[al,1]*self.atoms[0,at])*np.sqrt(self.mass[at])
                # ax-int coupling terms
                for gam in range(3):
                    Drot[3*at+gam,al] += np.sum([ I_gen_is[al,3+i]*np.sqrt(self.mass[at])*Bi_var[3*at+gam,i] for i in range(N_IR) ])
        
        ## internal rotations (eq. 12c)
        Dint = np.zeros((3*self.Nat,N_IR))
        for i in range(N_IR):
            for at in range(self.Nat):
                # ax-int coupling terms
                Dint[3*at+0,i] = (I_gen_is[3+i,1]*self.atoms[2,at]-I_gen_is[3+i,2]*self.atoms[1,at])*np.sqrt(self.mass[at])
                Dint[3*at+1,i] = (I_gen_is[3+i,2]*self.atoms[0,at]-I_gen_is[3+i,0]*self.atoms[2,at])*np.sqrt(self.mass[at])
                Dint[3*at+2,i] = (I_gen_is[3+i,0]*self.atoms[1,at]-I_gen_is[3+i,1]*self.atoms[0,at])*np.sqrt(self.mass[at])
                # int-int terms
                for ax in range(3):
                    Dint[3*at+ax,i] += np.sum([ I_gen_is[3+i,3+j]*np.sqrt(self.mass[at])*Bi_var[3*at+ax,j] for j in range(N_IR) ])
        
        ## I_gen is not diagonal, thus, Dtrans, Drot and Dint are coupling
        mixing = None #'block' # 'full'
        if mixing=='block':
            # block wise diagonalization of I_gen -> only if I_gen[ax,int]==0
            ovrlval, ovrlvec = np.linalg.eigh(I_gen[:3,:3]) # diagonalize overall subblock
            intval, intvec   = np.linalg.eigh(I_gen[3:,3:]) # diagonalize internal subblock
            # linear combination of Dtrans/Drot vectors to match the principle axes system
            for i in range(3):
                Dvar[:,0+i] = np.sum( [ovrlvec[j,i]*Dtrans[:,j] for j in range(3)] , axis=0)
                Dvar[:,3+i] = np.sum( [ovrlvec[j,i]*Drot[:,j] for j in range(3)] , axis=0)
            # linear combination of int rot vectors accoring diagonalization of I_gen int-int subblock
            # however, due to our choice of Bi_var, this block is already diagonal!
            for i in range(N_IR):
                Dvar[:,6+i] = np.sum( [intvec[j,i]*Dint[:,j] for j in range(N_IR)] , axis=0)
        
        elif mixing=='full':
            # full diagonalization of I_gen
            Igen_val, Igen_vec = np.linalg.eigh(I_gen) # is hermetian
            # !!!!!!!!
                
        else:
            # no mixing of the vectors
            Dvar[:,0:3] = np.copy(Dtrans)
            Dvar[:,3:6] = np.copy(Drot)
            Dvar[:,6:6+N_IR] = np.copy(Dint)
        
        #control = np.dot(Dvar.T,Dvar)[:6+N_IR,:6+N_IR]
        #print_mat(control)
        #plt.matshow(control, cmap='seismic',vmin=-np.abs(control.max()),vmax=np.abs(control.max()))
        #plt.matshow(control-np.identity(6+N_IR), cmap='seismic')
        #plt.colorbar()
        #plt.title('control matrix')
        #plt.show()
        
        ## for projector matrix (eq. 17)
        Proj = np.identity(3*self.Nat)
        for i in range(6+N_IR):
            Proj -= np.dot(Dvar[:,i:i+1],Dvar[:,i:i+1].T)
        # explicit formulation (eq. 16)
        #for i in range(3*Nat):
        #    for j in range(3*Nat):
        #        Proj[i,j] -= np.sum([Dvar[i,0+g1]*Dvar[j,0+g1] for g1 in range(3)])
        #        Proj[i,j] -= np.sum([Dvar[i,3+g2]*Dvar[j,3+g2] for g2 in range(3)])
        #        Proj[i,j] -= np.sum([Dvar[i,6+g3]*Dvar[j,6+g3] for g3 in range(N_IR)])
        #plt.matshow(Proj,cmap='seismic')
        #plt.colorbar()
        #plt.show()
        
        ## apply projector to Hessian to project out 6+N_IR motions (eq. 28)
        H_Mx = np.dot(self.sqMi, np.dot(self.Hess,self.sqMi)) / bohr2A**2
        H_Mx_proj = symm_mat(np.dot(Proj,np.dot(H_Mx,Proj))) # symmetrization irons out numerical errors
        
        ## diagonalize projected Hessian (section before eq. 30)
        Val_pr, Vec_Mpr = np.linalg.eigh(H_Mx_proj) # H_Mx_proj is hermetian, thus, use eigh
        sort = np.argsort(np.abs(Val_pr)) # in case, there are non-zero negative freqs
        Val_prs = Val_pr[sort] *eigval2s2
        Val_prs[:6+N_IR] = np.zeros(6+N_IR) # set first 6+N_IR eigvals to zero
        if any(Val_prs < 0.):
            Warning_flag = True
        else:
            Warning_flag = False
        Vec_Mprs = Vec_Mpr[:,sort]
        Dvar[:,6+N_IR:] = Vec_Mprs[:,6+N_IR:]
        w_pr = np.sign(Val_prs)*np.sqrt(np.abs(Val_prs))/(2.*np.pi*c_SI*100.) # in cm**-1
        
        Vec_pr = np.dot(self.sqMi,Vec_Mprs) # remove mass-weighting
        vib_redmass = np.zeros(self.DoV-N_IR)
        for i in range(self.DoV-N_IR):
            Ni = np.linalg.norm(Vec_pr[:,6+N_IR+i])
            vib_redmass[i] = Ni**(-2)
            Vec_pr[:,6+N_IR+i] /= Ni
        
        D_IR = np.dot(self.sqMi,Dvar[:,6:6+N_IR]) # remove mass-weighting
        IR_redmass = np.zeros(N_IR) # equal to eigvals of int-int subblock of I_gen
        for i in range(N_IR):
            Ni = np.linalg.norm(D_IR[:,i])
            IR_redmass[i] = (Ni)**(-2)
            D_IR[:,i] /= Ni
        
        ## compare Vec_pr with L and w_pr with vibs!!!
        # find matches between Vev_pr and L
        OL = np.abs( np.dot(Vec_pr[:,6+N_IR:].T,self.L).T )
        vib_matches = []
        for r in range(self.DoV):
            s_maxOL = np.argmax(OL[r,:])
            for s in range(self.DoV-N_IR):
                r_maxOL = np.argmax(OL[:,s])
                ol_max = OL[r_maxOL,s_maxOL]
                if ol_max == np.max(OL[r,:]) and ol_max == np.max(OL[:,s]):
                    vib_matches.append([r_maxOL,s_maxOL])
        
        vib_matches = np.array(vib_matches)
        nomatch = np.array([s for s in range(self.DoV-N_IR) if s not in vib_matches[:,0]])
        
        OL2 = np.zeros((len(nomatch),N_IR))
        IR_matches = []
        for i in range(len(nomatch)):
            for s in range(N_IR):
                OL2[i,s] = np.abs(np.dot(self.L[:,nomatch[i]],D_IR[:,s]))
            IR_matches.append([nomatch[i],np.argmax(OL2[i])])
        IR_matches = np.array(IR_matches)
        
        print('\nNormal coordinates with IntRots projected out')
        print('     full         projected        overlap')
        for i in range(self.DoV):
            pos = np.where(vib_matches[:,0] == i)[0]
            if len(pos) == 1:
                r,s = vib_matches[pos[0]]
                print('%2i: %7.2f     %2i: %7.2f         %5.3f' %(r,self.vibs_orig[r],s,w_pr[6+N_IR+s],OL[r,s] ))
            elif len(pos) == 0:
                pos2 = np.where(IR_matches[:,0] == i)[0]
                if len(pos2) == 1:
                    r,s = IR_matches[pos2[0]]
                    r2 = pos2[0]
                    print('%2i: %7.2f    IR%i: %7.2f         %5.3f' %(r,self.vibs_orig[r],s,w_IR[s],OL2[r2,s] ))
                elif len(pos2) == 0:
                    print('%2i: %7.2f     ' %(i, self.vibs_orig[i] ))
        
        
        ### safe normal coordinates
        f_s = open(self.directory+self.basen+'_NormCoords.txt', 'w')
        
        f_s.write(self.basen + ' info\n\n')
        f_s.write('N_IR %i\n\n' %N_IR)
        
        f_s.write('projected normal modes (mass-weighted):\n')
        for i in range(3*self.Nat):
            f_s.write('%+.6e  '*(self.DoV-N_IR) %(*Vec_Mprs[i,6+N_IR:],) + '\n' )
        f_s.write('\n')
        
        f_s.write('harmonic frequencies and corresponding red masses:\n')
        for i in range(self.DoV-N_IR):
            f_s.write('%10.5f    %8.5f' %(w_pr[6+N_IR+i],vib_redmass[i]) + '\n')
        f_s.write('\n')
        
        f_s.write('internal rotation coordinates (mass-weighted):\n')
        if N_IR > 0:
            for i in range(3*self.Nat):
                f_s.write('%+.6e  '*(N_IR) %(*Dvar[i,6:6+N_IR],) + '\n' )
        f_s.write('\n')
            
        f_s.write('Int Rot harmonic frequencies and corresponding red masses:\n')
        for i in range(N_IR):
            f_s.write('%10.5f    %8.5f' %(w_IR[i],IR_redmass[i]) + '\n')
        f_s.write('\n')
        
        f_s.close()
        
        np.savetxt(self.directory+self.basen+'_Projector.txt',Proj)
        np.savetxt(self.directory+self.basen+'_Dvar_mat.txt',Dvar)
        
        print('\nGeneralized inertia tensor:')
        print('    x         y         z     '+'    %i     '*N_IR %(*np.arange(N_IR),))
        print_mat(I_gen,fmt='%+9.3f')
        
        #
        if Warning_flag:
            print('\nWARNING! Projected Hessian has negative eigenvalues. Structure might not be a true minimum!')
            print('You can check the array w_pr for which are negative (the first 6+N_IR entries will be zero).')
            print('Plot the respective displacement using:        M.plot_NMD_3D(M.Vec_pr[:,6+M.N_IR+i])')
            print('where i is the index of the vibration (starting from 0).')
        
        ### plotting
        rot_ax_list = np.zeros((N_IR,3))
        ats_off_list = [0]*N_IR
        for i in range(N_IR):
            rot_ax_list[i], ats_off_list[i] = self.find_rotax(D_IR[:,i])
        
        # plot internal rotor displacement for inspection with bond_ats_save
        ssig = np.zeros(N_IR)
        for i in range(N_IR):
            if bond_ats_save[i] != 0.: # if there is a match in method1+2
                a1 = bond_ats_save[i][0]
                a2 = bond_ats_save[i][1]
                rot_ax_plot = np.array([self.atoms[:,a1],self.atoms[:,a2]]).T
                self.plot_NMD_3D(D_IR[:,i],rot_ax_plot,col='b')
            else:                                              # if there is no match   
                off = np.zeros(3) #atoms[:,ats_off_list[i]]
                rot_ax_plot = np.array([off,off+rot_ax_list[i]]).T
                self.plot_NMD_3D(D_IR[:,i],rot_ax_plot,col='b')
            
            # get sigma from user prompt
            ssig[i] = input('Enter rot. sigma of internal rotation %i: '%(i+1))
        
        return N_IR, w_pr[6+N_IR:], w_IR, Vec_pr, D_IR, ssig

    def find_rotax(self,D_IR):
        # find atoms with lowest displacements
        displ = np.zeros(self.Nat)
        for i in range(self.Nat):
            displ[i] = np.linalg.norm(D_IR[3*i:3*(i+1)])
        at1, at2 = np.argsort(displ)[:2] # lowest two displacements
        vec_i = self.atoms[:,at1] - self.atoms[:,at2] # initial guess for rotational axis
        vec_i /= np.linalg.norm(vec_i)
        
        # now minimize dot product with all displacements
        dot_sum = lambda x: np.abs(np.sum([ np.dot(x,D_IR[3*at:3*(at+1)]) for at in range(self.Nat) ]))
        cons = {'type': 'eq', 'fun': lambda x: np.linalg.norm(x)-1.} # constrain norm of x to 1
        res = minimize(dot_sum,vec_i, constraints=cons)
        
        return res.x, at1
    
    def IntRot_eval(self,ax_ats,at_list,sig):
        ''''ax_ats are two(!) atom indices defining the rotational axis,
        at_list are atom indices defining the rotating atoms,
        sig is the rotational symmetry number.'''
        # rotate in abc coords
        ovrlval, ovrlvec = np.linalg.eigh(self.II) # diagonalize overall subblock
        xyz_abc = np.dot(ovrlvec.T,self.atoms)
        
        ax = xyz_abc[:,ax_ats[1]]-xyz_abc[:,ax_ats[0]] # rotational axis
        ax_n = ax / np.linalg.norm(ax) # normalized rotational axis
        axoff = xyz_abc[:,ax_ats[0]] # starting point for rotational axis
        N_at = len(at_list) # number of rotating atoms
        
        # get coordinates and masses of IntRot
        xyz = np.array([ xyz_abc[:,at] for at in at_list ])
        masses = np.array([ self.mass[at] for at in at_list ])
        M_rot = np.sum(masses) # total mass of the rotor
        #print(xyz)
        #print(masses)
        #print(M_rot)
        
        # get CoM of rotor
        CoM = np.sum(masses*(xyz).T,axis=1)/M_rot
        
        # check, whether CoM is on rotational axis or not
        flag_asymm = False
        # calc projection of CoM onto rotational axis
        fac = np.dot(CoM-axoff, ax_n) # axoff + fac*ax_n gives projection of CoM onto rot axis
        projec = axoff + fac*ax_n
        # check distance between CoM and projection
        distance = np.linalg.norm(CoM - projec)
        if distance > 1e-2: # threshold whether CoM is ON rot axis or not
            flag_asymm = True
            print('Rotor is asymmetric! |CoM-ax_rot| = %.5f A' %distance)    
            print('CoM  is: %+.5f %+.5f %+.5f' %(*CoM,))
            print('Proj is: %+.5f %+.5f %+.5f' %(*projec,))
            #print('using x=%.4f' %fac)
        else:
            print('Rotor is symmetric! |CoM-ax_rot| = %.5f A' %distance) 
        
        # define axes of the rotor
        ez = np.copy(ax_n) # rotational axis is z-axis
        if flag_asymm: # x-axis goes through CoM of top
            P = np.copy(CoM)
        else: # x-axis is defined arbitrarily (since CoM lies ON z-axis)
            candidates = [x for x in at_list if x not in ax_ats] # atoms not on axis
            P = xyz_abc[:,candidates[-1]] # pick the last atom to define P
        fac = np.dot(P-axoff, ez) # axoff + fac*ez gives projection of P onto rot axis
        origin = axoff + fac*ez # origin of the coord system of the top
        pre_ex = P - origin # x-axis goes from origin to P (which is CoM or random atom)
        ex = pre_ex/np.linalg.norm(pre_ex)
        ey = np.cross(ez,ex)                       # or other way around?
        
        # check right-handedness
        # [Pitzer1946, p. 240] using the fact that the global coordinates are the canonical unit vectors
        alph_mat = np.array([ex,ey,ez])
        check = np.round(np.linalg.det(alph_mat),0)
        if check != 1.0:
            raise ValueError("Coordinate system of rotor is wrong! It's %.2f but should be +1" %check)
        
        # transform coordinates of rotor atoms into new axes
        xyz_rot = np.array([ np.dot(alph_mat,xyz[i]-origin) for i in range(N_at)])
        
        # define inertias in amu*A**2 [Pitzer1946, p. 240]
        A = np.sum([ masses[at]*(xyz_rot[at][0]**2 + xyz_rot[at][1]**2) for at in range(N_at) ])
        B = np.sum([ masses[at]*(xyz_rot[at][0] * xyz_rot[at][2]) for at in range(N_at) ])
        C = np.sum([ masses[at]*(xyz_rot[at][1] * xyz_rot[at][2]) for at in range(N_at) ])
        U = np.sum([ masses[at]*xyz_rot[at][0] for at in range(N_at) ])
        
        # beta vector [Pitzer1946, p. 240]
        alph = alph_mat.T
        beta = np.zeros(3)
        for i in range(3):
            A_term = alph[i,2]*A
            B_term = alph[i,0]*B
            C_term = alph[i,1]*C
            U_term = U*( alph[(i-1)%3,1]*origin[(i+1)%3] - alph[(i+1)%3,1]*origin[(i-1)%3] )
            beta[i] = A_term - B_term - C_term + U_term
        beta = beta
        
        # reduced moment of inertia
        Ival = np.copy(ovrlval) # overall moments of inertia
        Ir = A - np.sum( (alph[:,1]*U)**2/self.m_tot + beta**2/Ival )
        Ir_SI = Ir * amu2kg*1e-20
        
        # corresponding rotational constant in cm**-1
        Brot = h_SI/(8.*(np.pi)**2*c_SI*100.*Ir_SI)
        
        print('Ir  = %7.3f   amu*A**2' %Ir)
        print('B   = %9.5f cm**-1' %Brot)
        print('sig = %3.0f' %sig)
        
        return Ir, sig, Brot
    
    def distort(self,Li,scl):
        '''distort geometry along normal mode L (DoV array).
        scl is the step size in Angstrom (if Li is normalized!).'''
        xyz_dist = np.copy(self.atoms)
        for at in range(xyz_dist.shape[1]):
            xyz_dist[:,at] = self.atoms[:,at] + scl*Li[3*at:3*(at+1)]
        
        xyz_dist = xyz_dist.T
        fo=open(self.directory+self.basen+'_distort.xyz','w+')
        fo.write('%i\n' %self.Nat)
        fo.write(self.basen+' distorted\n')
        for i in range(self.Nat):
            fo.write(' %s  %.8f  %.8f  %.8f' %(self.atom_labs[i],*xyz_dist[i],) + '\n')
        fo.close()
        print('Distorted file written to %s%s_distorted.xyz' %(self.directory,self.basen))
        return
    
    def print_summary(self):
        print('\n\n= Summary of Molecule =\n')
        print('Number of atoms: %i' %self.Nat)
        print('Degrees of vibration:         %3i' %(self.DoV-self.N_IR))
        print('Degrees of internal rotation: %3i' %self.N_IR)
        print('Charge: %i, Multiplicity: %i' %(self.charge,self.multi))
        print('Eelec: %.8f Hartree' %self.Eelec)
        print('E_ZPE: %.8f Hartree' %self.ZPE)
        print('Mass: %.4f amu' %self.m_tot)
        if self.dipole != None:
            print('Dipole: %.4f Debye' %self.dipole)
        if self.polariz != None:
            print('Polarizability: %.4f Angstrom**3' %(self.polariz*1e30) )
        print('Rotational constants: %f %f %f GHz' %(self.RotA,self.RotB,self.RotC))
        print('Rotational symmetry number: %i' %self.sigma_OR)
        #print('Normal modes [cm-1]:')
        #print(*vibs)
        #print('\n')
    
    def write_mpr(self):
        '''Writes basen.mpr file. If Mobil=True the Q(l) data is also written into the .mpr file.'''
        print('\n== Writing %s.mpr file ==' %self.basen)
        
        fo = open(self.directory + self.basen +'.mpr', 'w+') # molecular properties
        fo.write('## data for ' + self.basen + '\n')
        fo.write('\n')
        fo.write('E_elec %.10f' %self.Eelec + '\n') # in Hartree
        fo.write('ZPE %.10f' %self.ZPE + '\n') # ZPE in Hartree
        fo.write('mass %.5f' %self.m_tot + '\n') # in amu
        fo.write('charge %.0f' %self.charge + '\n') # charge in elementary units
        fo.write('multiplicity %.0f' %self.multi + '\n') # multiplicity
        if self.dipole != None:
            fo.write('dipole_moment %.5f' %self.dipole + '\n') # multiplicity
        if self.polariz != None:
            fo.write('polarizability %.4e' %self.polariz + '\n') # multiplicity
        fo.write('\n')
        
        # format needs to be 'molec_{A_X_}a' where A is a string of a solvent
        # X is the number of that solvent and 'a' is a letter specifying the conformer
        try:
            Nstr = self.basen.split('_',1)[1].rsplit('_',1)[0]
            N = ''
            for Nstri in Nstr.split('_')[1::2]:
                N += Nstri + ' '
        except:
            N = '-1' # for solvents
        fo.write('Nclust ' + N + '\n')     # Nclust here!
        fo.write('\n')
        
        fo.write('RotConst/GHz ' + '%.8f  '*3 %(self.RotA,self.RotB,self.RotC) + '\n') # in GHz
        fo.write('RotSig %.0f' %self.sigma_OR + '\n')
        if self.dipole != None:
            fo.write('Rot_dipole %s' %self.dipole_ax + '\n')
        fo.write('\n')
        
        fo.write('IntRots (Ir/amu*A**2, w_harm/cm**-1, sigma)' + '\n')
        for i in range(self.N_IR):
            fo.write('%.6f    %.4f    %.0f' %(self.Iir[i],self.w_IR[i],self.sig_ir[i]) + '\n')
        fo.write('endIntRots' + '\n\n')
        
        fo.write('vibs/cm**-1' + '\n') # in cm^-1
        for i in range(self.DoV-self.N_IR):
            fo.write('%.4f' %self.w_pr[i] + '\n')
        fo.write('endvibs' + '\n')
        
        if hasattr(self, 'QDAT'):
            fo.write('\n')
            fo.write('# gst        q1st       q2st       q3st\n')
            for i in range(self.QDAT.shape[0]):
                fo.write('%.4e %.4e %.4e %.4e' %(*self.QDAT[i,:4],) + '\n')
        
        fo.close()
        print('Molecular Parameters written to ' + self.directory + self.basen +'.mpr')

#%%
direc = 'Systems/2-methyl-6-ethoxyquinoline_IPA/2_0_a/'
basen = '2-methyl-6-ethoxyquinoline_H2O_2_IPA_0_a'
#direc = 'Systems/solvents/IPA_2/'
#basen = 't-IPA'

M = molecule(direc,basen,EHigh=True,Mobil=True,IntRot=True)
M.print_summary()
M.write_mpr()

# move mpr file outside of local folder
os.replace(direc+basen+'.mpr', direc+'/../'+basen+'.mpr')

#%% look at vibrations after IntRots are projected out
# ii = 2 # which vibration (start counting at 0, up to M.DoV-M.N_IR)
# M.plot_NMD_3D(M.Vec_pr[:,6+M.N_IR+ii]*2.)
# print('w_pr = %.2f cm**-1' %(M.w_pr[ii]))