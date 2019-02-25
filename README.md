# MobCal-MPI
Parallelization of the commonly used MobCal suite to calculate ion mobilities and collision cross sections. From the Read_Me.txt: 

--------------------------------------------------- INTRODUCTION ---------------------------------------------------

This Python3 script is used to convert Gaussian .log outputs to MobCal .mfj inputs with MMFF94 atom labels and parameters. 
Steps 1 and 2 only need to be done the first time the script is used.

A PowerPoint file is included with this release that covers the fundamentals of the trajectory method, optimizations outlined in Analyst (2019), 144, 1660-1670, and file formats/ conversions used in the creation of MobCal inputs (.mjf).
While it is not necessary to read this, it is strongly recommended to understand what you are calculating, how the mfj_creator.py script works, and more importantly, how MobCal_MPI works.


------------------------------------------- INSTALLING RELEVANT PACKAGES -------------------------------------------

1. Installing OpenBabel

Install OpenBabel. This script has been optimized for use with OpenBabel v2.4.1, which can be found at http://openbabel.org/wiki/Main_Page
Corresponding reference: N M O'Boyle, M Banck, C A James, C Morley, T Vandermeersch, and G R Hutchison. "Open Babel: An open chemical toolbox." J. Cheminf. (2011), 3, 33.

2. Installing SDF2XYZ2SDF

Install sdf2xyz2sdf. This package can be obtained from:
Windows: https://sourceforge.net/projects/sdf2xyz2sdf/files/binaries/windows/  
Corresponding reference: P Tosco and T Balle. Journal of Molecular Modeling, 2011, 17, 3021-3023.

For simplicity later, please install this in the default directory specified by the sdf2xyz2sdf installer.

3. Add OpenBabel and sdf2tinkerxyz directories to PATH. 

If you specified the default install directories for both packages, these need to be added to PATH. The directories are: 
C:\Program Files\OpenBabel-2.4.1
C:\open3dtools\bin

---------------------------------------- GUASSIAN CALCULATION REQUIREMENTS ----------------------------------------

We recommend  that molecular structures input in to MobCal-MPI are optimized at the B3LYP/6-31++G(d,p) level of theory, using empirical corrections for dispersion interactions (empiricaldispersion=GD3) 
Partial charges are REQUIRED. For the highest accuracy, partial charges should be calculated using the ChelpG or MK partition scheme, and constrained to reproduce the molecular dipole moment. Utilization of Mulliken charges is associated with large errors in CCS calculation accuracy.
We recommend the pop=(mk,dipole) keyword if using Gaussian09/Gaussian16. 

For relevant discussion on the choice of partial charges, see Analyst (2019), 144, 1660-1670
An example input line for G09/G16 should look like: 

# opt freq b3lyp/6-31++g(d,p) int=ultrafine scf=xqc empiricaldispersion=GD3 pop=(mk,dipole)

------------------------------------ UTILIZING THE MFJ CREATOR PYTHON PACKAGE -------------------------------------

All functionality in the .mfj creator can be controlled in the mfj_creator.py file. Open the mfj_creator.py in IDLE by right-clicking the file, and selecting the option "Edit with IDLE"

1. Edit the "directory =" field (line 3) to the directory which contains your Gaussian .log files. Ensure the path is wrapped by r'/filedirectory/here'. An example input line looks like:
directory = r'C:\Users\JamFlex\ModJed\Mobcal_ready'

2. There is an option to convert all .log files in your directory to MobCal inputs (.mfj). To do so, change the "csv_list =" (line 8) variable to:
csv_list = r''

If you wish to only convert select .log files in the directory to MobCal inputs (.mfj), you may do so by creating a .csv that contains the names of the .log files you wish to convert.
To specify a .csv, add the path and file name of the .csv to line 8. Ensure the path is wrapped by r'/filedirectory/csv_file.csv'. An example input line looks like:
csv_list = r'C:\Users\JamFlex\ModJed\Mobcal_ready\files.csv'

3. If sdf2xyz2sdf was not installed in the default directory, you can specify the location of it on line 12. If you installed in the default directory, you should see it in the following directory:

sdf2xyz2sdf_Directory = r'C:\open3dtools\bin\sdf2tinkerxyz.exe'

Again only change the sdf2xyz2sdf_Directory if you did not install in the default directory(C:\open3dtools\bin\sdf2tinkerxyz.exe)
DO NOT remove the r' before the directory or the ' trailing the directory!

4. You can define partial charge schemes on line 15 using the "charge =" variable. There are 3 options:

charge = 'calc' 
This specifies that a unique partial charge scheme is going to be read. Currently this MFJ creator can read charges predicted by electrostatic potential mapping methods (eg. ChelpG, MK) by Gaussian directly from the .log output file. 
Recommended. 

charge = 'equal' 
This specifies that equivalent partial charges are assigned to each atom by the relationship partial charge = (ion charge) / (number of atoms)
Not recommended. 

charge = 'none' 
Assign all atomic charges as zero, which effectively removes the ion-induced dipole and ion-quadrupole terms for the ion-gas potential. 
Not recommended. 

********
Note that regardless of charge schemes, this script was built in such a way such that Gaussian .log files must contain partial charges calculated by either the MK/ChelpG schemes. 
In addition, in order for the .mfj to be read correctly by MobCal-MPI, it much contain partial charges (even if a charge scheme of equal or none is specified). Without these, the code will not function (this is intentional!)
********

5. You can define parameters to use in the trajectory method calculations which dictate how many trajectories are performed. 
Respective nomenclature is number of cycles (itn), number of velocity point integrations (inp), and number of impact parameter integrations (imp). 
The total number of trajectories run will be itn*inp*imp. 

The number of cycles (itn) is defined on line 18. The CCS reported is taken as the average across the number of cycles, for which the recommended value is 10. Your input line should look like:

cycles = 10

The number of velocity point integrations (inp) is defined on line 21. The recommended value is 48. Your input should look like:

v_integrations = 48

The number of impact parameter integrations (imp) is defined on line 24. The recommended value is 512. Your input should look like:

b_integrations = 512

********
Note that the number of velocity and impact parameter integration points must scale accordingly to the number of core you intend to run this method on. Thus, it is required that you specify an integer number of cores that can evenly distribute trajectory calculations across all cores. 
While this does not affect input file creation, it is meant to be an error check. Typically, Mobcal-MPI performs most efficiently on 8 or 16 cores. Thus, imp and inp must be divisible by the number of cores. 

The number of cores is specified on line 27 and must be an integer. Your input should look like: 

n_cores = 16
********

6. You can define the gas you wish to run trajectory method CCS calculations in. This is done on line 30. Currently, support is only available for He and N2. 

For He, your input should look like: 

gas = 'He'

For N2, your input should look like: 

gas = 'N2'

7. Please do not edit any other fields or portions of the code. 

8. Execute the code by pressing F5 (if you are using IDLE). You may be prompted to save the mfj_creator.py. Hit okay. 

After mfj creation is complete, we recommend that manual verification of atom type assignments is completed. 
If an error message is received where an unidentified atom type has been identified, this must be input manually. This can be done by correlating atom type in the atom_type.prm file, finding the atom class, and correlating that with the vdw parameters found in vdw.prm. 
For a description of these files, please see the included PowerPoint. 

------------------------------------------ THE MOBCAL INPUT FILE (.mfj) ------------------------------------------

Upon successful conversion, your .mfj input will have the following format:

Title
Number of conformations (always 1; for multiple conformers, we recommend creation of separate .mfj inputs)
Number of atoms
Units for xyz coordinates (ang for angstroms, au for atomic units)
Charge scheme (calc for non-uniform, custom charges; equal for equal charges (i.e., 1 / N(atoms)); none for no partial charges)
Scaling factor for xyz coordinates (typically always 1.0000)
itn inp imp gas# (1 for He, 2 for N2)
x	y	z	mass	partial charge	alpha_i	Ai	Ni	Gi

Example:

AMIFOSTINE
1
28
ang
calc
1.0000
10 48 512 2
  1.290700	     0.833300	    -1.273600	    31.972	   -0.193666	   3.00	 4.800	3.320	   1.345
  1.835900	    -0.567900	     0.215300	    30.974	    0.928334	  1.600	 4.500	3.320	   1.345
  1.601500	     0.228400	     1.605900	    15.995	   -0.608579	   0.70	 3.150	3.890	   1.282
  3.441700	    -0.700200	     0.177400	    15.995	   -0.621605	   0.70	 3.150	3.890	   1.282
  1.123800	    -1.872800	     0.078400	    15.995	   -0.623803	   0.75	 3.150	3.890	   1.282
 -1.323000	     0.767800	     0.692000	    14.003	   -0.039906	   1.15	 2.820	3.890	   1.282
 -1.458700	    -1.483100	    -0.811300	    14.003	   -0.302638	   1.00	 2.820	3.890	   1.282
 -2.757700	     0.544900	     1.014800	    12.000	   -0.136973	  1.050	 2.490	3.890	   1.282
 -3.480200	    -0.274000	    -0.071900	    12.000	   -0.019077	  1.050	 2.490	3.890	   1.282
 -1.054200	     1.971400	    -0.111000	    12.000	   -0.040526	  1.050	 2.490	3.890	   1.282
  0.439900	     2.194000	    -0.338400	    12.000	   -0.014244	  1.050	 2.490	3.890	   1.282
 -2.860900	    -1.659300	    -0.288600	    12.000	    0.095374	  1.050	 2.490	3.890	   1.282
 -2.799000	    -0.005900	     1.960400	     1.008	    0.091435	  0.250	 0.800	4.200	   1.209
 -3.284100	     1.496100	     1.171900	     1.008	    0.080214	  0.250	 0.800	4.200	   1.209
 -3.500200	     0.275800	    -1.022000	     1.008	    0.046225	  0.250	 0.800	4.200	   1.209
 -4.525400	    -0.406500	     0.225100	     1.008	    0.065223	  0.250	 0.800	4.200	   1.209
 -1.458100	     2.875500	     0.371900	     1.008	    0.042858	  0.250	 0.800	4.200	   1.209
 -1.565300	     1.866600	    -1.074800	     1.008	    0.082288	  0.250	 0.800	4.200	   1.209
  0.974400	     2.358200	     0.598800	     1.008	    0.076472	  0.250	 0.800	4.200	   1.209
  0.589100	     3.075000	    -0.968200	     1.008	    0.131005	  0.250	 0.800	4.200	   1.209
 -3.432300	    -2.270100	    -0.990200	     1.008	    0.053687	  0.250	 0.800	4.200	   1.209
 -2.779200	    -2.204400	     0.655600	     1.008	    0.072372	  0.250	 0.800	4.200	   1.209
 -1.429400	    -1.306500	    -1.816500	     1.008	    0.278131	  0.150	 0.800	4.200	   1.209
 -0.795400	    -2.237000	    -0.596000	     1.008	    0.326918	  0.150	 0.800	4.200	   1.209
  2.343200	     0.155300	     2.228600	     1.008	    0.491602	  0.150	 0.800	4.200	   1.209
  3.745900	    -1.567900	    -0.135000	     1.008	    0.495074	  0.150	 0.800	4.200	   1.209
 -0.776100	     0.804000	     1.549800	     1.008	    0.178081	  0.150	 0.800	4.200	   1.209
 -1.083700	    -0.620300	    -0.299900	     1.008	    0.065725	  0.150	 0.800	4.200	   1.209

--------------------- COMPILING MOBCAL MPI AND PERFORMING TRAJECTORY METHOD CCS CALCULATIONS ---------------------

MobCal-MPI was developed for use on high performance computing architectures, many of which are equipped with legacy Fortran compilers. 

To compile MobCal-MPI:

1. Place the MobCal_MPI.f code in a directory on the parallel platform. 

2. Compile using the following command (or using an alternative legacy Fortran compiler): 
mpif90 -o MobCal_MPI.exe MobCal_MPI.f -Ofast 

3. A bash script for SLURM schedulers has been provided (runmob_MPI.sh). To use this:

Place script in the same directory containing your .mjf files and compiled MobCal_MPI (will be a .exe extension after compilation)

Change file permission to allow for execution using the following command: 
chmod a+x runmob_MPI.sh

Execute the bash script to submit all .mfj files in that directory to MobCal_MPI using the following command:
./runmob_MPI.sh

Use this for submission on all SLURM schedulers. 

------------------------------------------------ GET IN TOUCH ---------------------------------------------------

Questions? Comments? Concerns? Email us! 

cieritano@uwaterloo.ca
shopkins@uwaterloo.ca


