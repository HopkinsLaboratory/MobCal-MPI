#!/bin/bash

# set ion name:
ION_NAME="Acetone"
ION_CID="180"

# set MobCal parameters
T_MOBCAL="300 1000 81"
ITN_MOBCAL="10"
INP_MOBCAL="192"
IMP_MOBCAL="1024"

# set main folder
MAIN_DIREC="/home/$USER/"

######## end of user input ########

# download chemical and rename
make download-chemical CHEMICAL_CID=$ION_CID
cd "${MAIN_DIREC}Results/"
mv "COMPOUND_CID_${ION_CID}" "${ION_NAME}"
cd "${ION_NAME}/Chemical_Structure"
mv *"_${ION_CID}.sdf" "${ION_NAME}.sdf"

# create output folders
cd ..
mkdir -p "Ions/"
mkdir -p "Conformers/"
mkdir -p "ORCA/Inputs/"
mkdir -p "ORCA/Outputs/"
mkdir -p "MobCal/Inputs/"
mkdir -p "MobCal/Outputs/"

# create scr folders
cd "${MAIN_DIREC}MobCal-MPI_21/src/build"
mkdir -p "$ION_NAME"
cd "$ION_NAME"
mkdir -p ions
mkdir -p positive_build
mkdir -p orca_outputs
mkdir -p mobcal_inputs
mkdir -p mobcal_outputs

# exit to MobCal folder for make execution
cd "${MAIN_DIREC}MobCal-MPI_21"

# define directory variables
CUR_MAKEFILE="${MAIN_DIREC}MobCal-MPI_21/Makefile"
INPUT_PATH="${MAIN_DIREC}Results/${ION_NAME}/Chemical_Structure/${ION_NAME}.sdf"

ION_BUILD_DIRECTORY="${MAIN_DIREC}MobCal-MPI_21/src/build/${ION_NAME}/"
POS_ION_BUILD_DIRECTORY="${ION_BUILD_DIRECTORY}positive_build/"

ION_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/Ions/"
CONFORMER_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/Conformers/"
ORCA_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/ORCA/"
MOBCAL_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/MobCal/"

# CREST: create ions and sample conformers of monomer and dimer
mkdir -p "${POS_ION_BUILD_DIRECTORY}conjugate_build"
make --ignore-errors --makefile=$CUR_MAKEFILE --directory=$POS_ION_BUILD_DIRECTORY generate-conjugate INPUT_FILE=$INPUT_PATH OUTPUT_DIRECTORY=$ION_OUTPUT_DIRECTORY CHARGE="1"

CONFORMER_DIRECTORY="${ION_BUILD_DIRECTORY}ions/${ION_NAME}_protonate.xyz/"
mkdir -p $CONFORMER_DIRECTORY
make --makefile=$CUR_MAKEFILE --directory=$CONFORMER_DIRECTORY generate-conformers SIMILARITY_THRESHOLD=95 INPUT_FILE="${ION_OUTPUT_DIRECTORY}${ION_NAME}_protonate.xyz" OUTPUT_DIRECTORY=$CONFORMER_OUTPUT_DIRECTORY CHARGE="1" MULTIPLICITY="1"

# ORCA: create inp files and run them
make --makefile=$CUR_MAKEFILE --directory="${CONFORMER_DIRECTORY}split/uniques" generate-orca-inputs CHARGE="1" MULTIPLICITY="1" INPUT_DIRECTORY="." OUTPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Inputs/"

make --makefile=$CUR_MAKEFILE --directory="${ION_BUILD_DIRECTORY}orca_outputs/" generate-orca-outputs INPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Inputs/" OUTPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Outputs/"

# MobCal: create mfj files and run them
make --makefile=$CUR_MAKEFILE --directory="${ION_BUILD_DIRECTORY}mobcal_inputs/" generate-mobcal-inputs INPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Outputs/" OUTPUT_DIRECTORY="${MOBCAL_OUTPUT_DIRECTORY}Inputs/" CHARGE_SETTING="calc" MFJ_CYCLES=$ITN_MOBCAL MFJ_VELOCITY_INTEGRATION=$INP_MOBCAL MFJ_IMPACT_INTEGRATION=$IMP_MOBCAL MFJ_CARRIER_GAS=2 MFJ_TEMPERATURE="${T_MOBCAL}" MFJ_EMPIRICAL_CORRECTION=0

make --makefile=$CUR_MAKEFILE --directory="${ION_BUILD_DIRECTORY}mobcal_outputs/" generate-mobcal-outputs INPUT_DIRECTORY="${MOBCAL_OUTPUT_DIRECTORY}Inputs/" OUTPUT_DIRECTORY="${MOBCAL_OUTPUT_DIRECTORY}Outputs/"

# cleanup
make --makefile=$CUR_MAKEFILE --directory="${MAIN_DIREC}MobCal-MPI_21/src/build/" remove-build-directory DIRECTORY=$ION_BUILD_DIRECTORY
