# MobCal-MPI

MobCal-MPI is a parallelized version of the widely-used MobCal suite developed by the [Jarrold group](https://pubs.acs.org/doi/10.1021/jp961623v), further modified to use atom-specific interaction potentials and two-temperature theory for efficient and accurate calculation of ion mobilities and collision cross sections. The 'X' module (MobCal-MPI-X) integrates with [CREST](https://crest-lab.github.io/crest-docs/) and [ORCA](https://orcaforum.kofo.mpg.de) quantum-chemical calculation packages to automate the CCS calculation process. 

MobCal-MPI is free for academic use. We only ask that you cite its corresponding publications:

[Analyst 2023, 148 (14), 3257–3273](https://doi.org/10.1039/D3AN00545C).

[Analyst, 2019,144, 1660-1670](https://doi.org/10.1039/C8AN02150C). 

[Source code and user manual)](https://zenodo.org/records/11426097)

## Fortran code Compilation Instructions

Detailed compilation instructions on high-performance computing (HPC) and standalone Unix platforms are provided in the [MobCal-MPI manual](/Manual/MobCal-MPI_User_manual.pdf). For the convenience of users who do not have access to HPC architecture, we have also provided a [Quick Guide](/Manual/Quick%20Guide%20for%20installing%20Intel%20OneAPI%20to%20compile%20Fortran%20code.pdf) for installing Intel OneAPI and the associate mpiifort compiler to local Unix systems.

## Automated CCS calculations (MobCal-MPI-X)
Details to follow. Stay tuned.

## GUI 
A [GUI](/GUI/Launcher.py) is provided to streamline the creation and analysis of MobCal-MPI files. Detailed documentation is provided in the [manual](/Manual/MobCal-MPI_User_manual.pdf), although for users who desire a QuickStart workflow, the GUI requires the following prerequisites to load:

- **Python 3.12+**: Download and install from [Python Downloads](https://www.python.org/downloads/).
- **GitHub Desktop**: Download from [GitHub Desktop](https://desktop.github.com/). Please install Git to the default location and add to PATH if prompted to do so. A GitHub account is not required for installation.
- **Git**: Download from [Git SCM](https://git-scm.com/). Please install Git to the default location and add to PATH if prompted to do so. A GitHub account is not required for installation.
- **OpenBabel 2.4.1**: Download and install from [sourceforge](https://sourceforge.net/projects/openbabel/files/openbabel/2.4.1/) to the default location and add the download location to your system's PATH.
- **SDF2XYZ2SDF**: Download and install from [sourceforge](https://sourceforge.net/projects/sdf2xyz2sdf/files/binaries/windows/) to the default location and add the download location to your system's PATH.

The GUI also requires several Python packages. Install them using the following command in your command prompt:

```console
pip install numpy pyqt6 scipy matplotlib gitpython
```

Once installed, the GUI can be run using your preferred Python environment. 

## Contact Information
Questions, comments, or concerns? Feel free to reach out via email or GitHub Issues:

- [a2haack@uwaterloo.ca](mailto:a2haack@uwaterloo.ca)
- [cieritano@uwaterloo.ca](mailto:cieritano@uwaterloo.ca)

