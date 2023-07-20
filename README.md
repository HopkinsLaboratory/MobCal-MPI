# MobCal-MPI
Parallelization of the commonly used MobCal suite to calculate ion mobilities and collision cross sections (v1.2). For full documentation, see the MobCal-MPI user guide.

# Changes in v2.0 <h2>
The latest release of MobCal-MPI (v2.0) expands on its predecessor by implementing two-temperature theory to accuately compute ion mobilities and CCSs at arbitrary field strengths. All functionality from v1.2 is retained, including calulcation of CCSs within the low-field limit.

A summary of the changes are included below:

-	Implementation of two-temperature theory (2TT) within the MobCal-MPI framework 
-	Implementation of an empirical correction to 2TT, thus correcting systematic deviations between calculated and experimental high-field mobilities 
-	Updated the methodology for evaluating uncertainty in the CCS calculated from the trajectory method
-	Changed the nature of the grid in which ion-neutral relative velocities are sampled during the simulation of collision dynamics from a weighted spacing to a linear spacing. 
-	Fixed a bug where the temperature in which N2 rotational orientations were weighted was fixed at 500 K rather than set to the user-defined temperature. 

For a full description of each change, we direct users to the manuscript accompanying MobCal-MPI 2.0. 

Analyst 2023, 148 (14), 3257–3273. https://doi.org/10.1039/D3AN00545C.

# Changes in v1.2 <h2>

This release of MobCal-MPI (v 1.2) is equipped with a Python-based Graphical User Interface (GUI) for:


1.Conversion of Gaussian output files (.log) into MobCal-MPI input files (.mfj)

2.Calculation of Boltzmann weights of isomer populations for the generation of Boltzmann-weighted CCSs

3.Extraction of CCSs from MobCal-MPI output files (.mout)

4.Multiple temperatures can now be specified in a single .mfj  input file

5.The reporting of CCS now occurs at the end of each .mout file in a formatted table 



**Get in touch**

Questions? Comments? Concerns? Email us! 

a2haack@uwaterloo.ca
cieritano@uwaterloo.ca
shopkins@uwaterloo.ca



