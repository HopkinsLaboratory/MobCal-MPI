# MobCal-MPI 2.0.1

MobCal-MPI is a parallelized version of the widely-used MobCal suite, designed for calculating ion mobilities and collision cross sections (2.0).

## Compilation Instructions

Detailed compilation instructions on high-performance computing (HPC) and standalone Unix platforms are provided in the [MobCal-MPI 2.0.1 manual](/Manual/MobCal-MPI_User_manual.pdf). For the convenience of users who do not have access to HPC architecture, we have also provided a [Quick Guide](/Manual/Quick%20Guide%20for%20installing%20Intel%20OneAPI%20to%20compile%20Fortran%20code.pdf) for installing Intel OneAPI and the associate mpiifort compiler to local Unix systems.

## Changes in v2.0.1
The latest release, v2.0.1, builds upon its predecessor by incorporating two-temperature theory (2TT) to accurately compute ion mobilities and CCSs at arbitrary field strengths. It addresses fixes to the GUI and some compilation errors. The release also includes a detailed manual covering the GUI's functionality, compiling the MobCal-MPI Fortran code on various Unix platforms, and an Appendix providing a step-by-step tutorial on calculating an ion's CCS from scratch.

### Summary of Changes:
- Implementation of two-temperature theory (2TT) within the MobCal-MPI framework
- Empirical correction to 2TT to rectify systematic deviations between calculated and experimental high-field mobilities
- Updated methodology for evaluating uncertainty in CCS calculated from the trajectory method
- Changed the nature of the grid for ion-neutral relative velocities sampling from weighted spacing to linear spacing
- Fixed a bug where the temperature for N2 rotational orientations weighting was fixed at 500 K instead of being user-defined

For a comprehensive description of each change, refer to the accompanying manuscript: [Analyst 2023, 148 (14), 3257–3273](https://doi.org/10.1039/D3AN00545C).

#### Changes from v1.2
This release (v1.2) of MobCal-MPI is equipped with a Python-based Graphical User Interface (GUI) for:

- Conversion of Gaussian output files (.log) into MobCal-MPI input files (.mfj)
- Calculation of Boltzmann weights of isomer populations for the generation of Boltzmann-weighted CCSs
- Extraction of CCSs from MobCal-MPI output files (.mout)
- Multiple temperatures can now be specified in a single .mfj input file
- Reporting of CCS now occurs at the end of each .mout file in a formatted table

## Contact Information
Questions, comments, or concerns? Feel free to reach out via email or GitHub Issues:

- [a2haack@uwaterloo.ca](mailto:a2haack@uwaterloo.ca)
- [cieritano@uwaterloo.ca](mailto:cieritano@uwaterloo.ca)
- [shopkins@uwaterloo.ca](mailto:shopkins@uwaterloo.ca)

