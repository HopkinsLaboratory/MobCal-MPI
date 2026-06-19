#!/bin/bash

# define list of CIDs and Names
numberToName=(
  [8461]="24DNT"
  [8376]="TNT"
  [6954]="PicricAcid"
  [10178]="Tetryl"
  [6518]="PETN"
  [5284553]="ETN"
)

# set MobCal parameters
T_MOBCAL="353 800 81"
ITN_MOBCAL="10"
INP_MOBCAL="192"
IMP_MOBCAL="1024"

# set main folder
MAIN_DIREC="/home/$USER/"

######## end of user input ########

# iterate through list:
for ION_CID in "${!numberToName[@]}"; do
  # Get the in name corresponding to the current CID
  ION_NAME=${numberToName[$ION_CID]}
  
  # echo
  echo " "
  echo "==> Processing $ION_NAME (CID $ION_CID) <=="
  
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
  mkdir -p negative_build
  mkdir -p orca_outputs
  mkdir -p mobcal_inputs
  mkdir -p mobcal_outputs
  
  # exit to MobCal folder for make execution
  cd "${MAIN_DIREC}MobCal-MPI_21"
  
  # define directory variables
  CUR_MAKEFILE="${MAIN_DIREC}MobCal-MPI_21/Makefile"
  INPUT_PATH="${MAIN_DIREC}Results/${ION_NAME}/Chemical_Structure/${ION_NAME}.sdf"
  
  ION_BUILD_DIRECTORY="${MAIN_DIREC}MobCal-MPI_21/src/build/${ION_NAME}/"
  NEG_ION_BUILD_DIRECTORY="${ION_BUILD_DIRECTORY}negative_build/"
  
  ION_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/Ions/"
  CONFORMER_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/Conformers/"
  ORCA_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/ORCA/"
  MOBCAL_OUTPUT_DIRECTORY="${MAIN_DIREC}Results/${ION_NAME}/MobCal/"
  
  # CREST: create radical anion and sample its conformers
  #/usr/local/bin/obabel -isdf $INPUT_PATH -oxyz --gen3D > "${NEG_ION_BUILD_DIRECTORY}${ION_NAME}.xyz"
  mkdir -p "${NEG_ION_BUILD_DIRECTORY}ion_build"
  make --ignore-errors --makefile=$CUR_MAKEFILE --directory=$NEG_ION_BUILD_DIRECTORY generate-ion INPUT_FILE=$INPUT_PATH OUTPUT_DIRECTORY=$ION_OUTPUT_DIRECTORY CHARGE="-1"
  
  mv "${ION_OUTPUT_DIRECTORY}${ION_NAME}.sdf" "${ION_OUTPUT_DIRECTORY}${ION_NAME}_radical-anion.xyz"
  
  CONFORMER_DIRECTORY="${ION_BUILD_DIRECTORY}ions/${ION_NAME}_radical-anion.xyz/"
  mkdir -p $CONFORMER_DIRECTORY
  make --makefile=$CUR_MAKEFILE --directory=$CONFORMER_DIRECTORY generate-conformers SIMILARITY_THRESHOLD=95 INPUT_FILE="${ION_OUTPUT_DIRECTORY}${ION_NAME}_radical-anion.xyz" OUTPUT_DIRECTORY=$CONFORMER_OUTPUT_DIRECTORY CHARGE="-1" MULTIPLICITY="2"
  
  # ORCA: create inp files and run them
  make --makefile=$CUR_MAKEFILE --directory="${CONFORMER_DIRECTORY}split/uniques" generate-orca-inputs CHARGE="-1" MULTIPLICITY="2" INPUT_DIRECTORY="." OUTPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Inputs/"
  
  make --makefile=$CUR_MAKEFILE --directory="${ION_BUILD_DIRECTORY}orca_outputs/" generate-orca-outputs INPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Inputs/" OUTPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Outputs/"
  
  # MobCal: create mfj files and run them
  make --makefile=$CUR_MAKEFILE --directory="${ION_BUILD_DIRECTORY}mobcal_inputs/" generate-mobcal-inputs INPUT_DIRECTORY="${ORCA_OUTPUT_DIRECTORY}Outputs/" OUTPUT_DIRECTORY="${MOBCAL_OUTPUT_DIRECTORY}Inputs/" CHARGE_SETTING="calc" MFJ_CYCLES=$ITN_MOBCAL MFJ_VELOCITY_INTEGRATION=$INP_MOBCAL MFJ_IMPACT_INTEGRATION=$IMP_MOBCAL MFJ_CARRIER_GAS=2 MFJ_TEMPERATURE="${T_MOBCAL}" MFJ_EMPIRICAL_CORRECTION=0
  
  make --makefile=$CUR_MAKEFILE --directory="${ION_BUILD_DIRECTORY}mobcal_outputs/" generate-mobcal-outputs INPUT_DIRECTORY="${MOBCAL_OUTPUT_DIRECTORY}Inputs/" OUTPUT_DIRECTORY="${MOBCAL_OUTPUT_DIRECTORY}Outputs/"
  
  # cleanup
  make --makefile=$CUR_MAKEFILE --directory="${MAIN_DIREC}MobCal-MPI_21/src/build/" remove-build-directory DIRECTORY=$ION_BUILD_DIRECTORY
  
done