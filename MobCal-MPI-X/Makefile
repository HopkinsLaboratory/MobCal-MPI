# Set shell to bash
SHELL:=bash

# Pass down the makefile list
export MAKEFILE_LIST
CUR_MAKEFILE:=$(abspath $(lastword $(MAKEFILE_LIST)))

# Pass down build directory cleanup flag
CLEANUP_BUILD=""
export CLEANUP_BUILD

# Determine the appropriate calling location of the Makefile
CALL_DIRECTORY:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

# Generate file paths to the scripts and chemical structure directories
SCRIPTS_DIRECTORY:=$(CALL_DIRECTORY)/src/
BUILD_DIRECTORY:=$(SCRIPTS_DIRECTORY)build/
FORTRAN_DIRECTORY:=$(CALL_DIRECTORY)/
RESULTS_DIRECTORY:=$(abspath $(CALL_DIRECTORY)/../Results)/
RESULTS_CHEMICAL_STRUCTURE_DIRECTORY:=Chemical_Structure/
RESULTS_IONS_OUTPUT_DIRECTORY:=Ions/
RESULTS_CONFORMERS_OUTPUT_DIRECTORY:=Conformers/
RESULTS_ORCA_INPUTS_OUTPUT_DIRECTORY:=ORCA/Inputs/
RESULTS_ORCA_OUTPUTS_OUTPUT_DIRECTORY:=ORCA/Outputs/
RESULTS_MOBCAL_INPUTS_OUTPUT_DIRECTORY:=MobCal/Inputs/
RESULTS_MOBCAL_OUTPUTS_OUTPUT_DIRECTORY:=MobCal/Outputs/

# Variables to override for SMILES input
CHEMICAL_STRING:=""
CHEMICAL_NAME:=""
CHEMICAL_CID:=""

# Variables to override for ion generation
INPUT_FILE:=

# Variable to override for input check
INPUT_VALUE:=
INPUT_VARIABLE:="Not provided"

# Variables to override for PubChem data processing
PUBCHEM_DATA:=
OUTPUT_DIRECTORY:=
PUBCHEM_DATA_OUTPUT_FILE:="README.md"

# Common programs/scripts
MPIRUN:=mpirun
FORTRAN_COMPILER:=mpifort
OPEN_BABEL:=$(shell which obabel)
ORCA:=$(shell which orca)
SDF2TINKERXYZ:=$(shell which sdf2tinkerxyz)
EMPIRICAL_CORRECTIONS_CHECK=$(SCRIPTS_DIRECTORY)determine_empirical_correction.sh
ION_SPLITTER:=python3 $(SCRIPTS_DIRECTORY)ion-library-file-splitter.py
COMPARE_MASSES:=python3 $(SCRIPTS_DIRECTORY)compare_masses.py
GJF_ORCA_CONVERSION:=python3 $(SCRIPTS_DIRECTORY)gjf_Orca_conversion.py
THERMO_CHEM_EXTRACTOR:=python3 $(SCRIPTS_DIRECTORY)orca_thermochem_extractor.py
MFJ_CREATOR:=python3 $(SCRIPTS_DIRECTORY)mfj_creator.py
MOBCAL_OUTPUT_ANALYZER:=python3 $(SCRIPTS_DIRECTORY)GUI/output-analyzer/mout_Analyzer.py
DIMER_APPEND:=python3 $(SCRIPTS_DIRECTORY)append_solvent_chemical.py
CORE_COUNTER=$(SCRIPTS_DIRECTORY)/avail_core.sh
MEM_COUNTER=$(SCRIPTS_DIRECTORY)/avail_mem.sh
XYZ_FILE_FORMAT_CHECKER=python3 $(SCRIPTS_DIRECTORY)/check_xyz_file.py
CHECK_CHARGE_STRING=$(SCRIPTS_DIRECTORY)/check_string.sh
CHECK_MULTIPLICITY_STRING=$(SCRIPTS_DIRECTORY)/check_multiplicity.sh
CHECK_COSINE_SIMILAIRITY=$(SCRIPTS_DIRECTORY)/check_cosine.sh
CURRENT_TIMESTAMP_COMMAND=date +%s

# Variables to hold git branches
CURRENT_BRANCH:=$(shell git branch | grep "*" | sed "s/* //")
MAIN_BRANCH:=master
UPDATE_MESSAGE:="Automated ion library commit message."
UPDATE_BRANCH:=$(shell date | sed "s/ /_/g" | sed "s/:/_/g")
MAKEFILE_LOG:=makefile.illog

# Variables to hold values for chemical processing
DIVIDER_COUNT=$(shell printf "%s\n" $(shell stty size) | tail -n 1)
AVAILABLE_CORES:=$(shell $(CORE_COUNTER))
AVAILABLE_RAM:=$(shell $(MEM_COUNTER))
CONFORMERS_FILE:=crest_conformers.xyz
ION_SPLITTER_OUTPUT_DIRECTORY:=split
MASS_COMPARISON_OUTPUT_DIRECTORY:=uniques
MFJ_CREATOR_OUTPUTS_DIRECTORY=Mobcal_Inputs
FORTRAN_FLAGS:=-Og
MOBCAL_INPUT_FILE:=mobcal_input.run
DIMER_FILE=dimer.xyz
MFJ_CREATOR_LOG=mobcal.mfjcreatorlog
PIPELINE_SIMILARITY_THRESHOLD=95
PIPELINE_MFJ_CHARGE_SETTING="calc"
PIPELINE_MFJ_CYCLES=10
PIPELINE_MFJ_VELOCITY_INTEGRATION=200
PIPELINE_MFJ_IMPACT_INTEGRATION=1024
PIPELINE_MFJ_CARRIER_GAS=2
PIPELINE_MFJ_TEMPERATURE="298 1008 71"

update:
	@git stash &>> $(MAKEFILE_LOG)
	@git checkout $(MAIN_BRANCH) &>> $(MAKEFILE_LOG)
	@git pull --ff origin $(MAIN_BRANCH) &>> $(MAKEFILE_LOG)
	@git checkout $(CURRENT_BRANCH) &>> $(MAKEFILE_LOG)
	@STASH_RESPONSE=`git stash show 2> /dev/null`; if ! [[ -z "$$STASH_RESPONSE" ]]; then git stash pop &>> $(MAKEFILE_LOG); fi

process-chemical:
	@make --no-print-directory --makefile=$(CUR_MAKEFILE) check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	$(eval INPUT_FILE_BASENAME=$(basename $(notdir $(shell printf "%s" "$(INPUT_FILE)" | sed "s/ /_/g"))))
	$(eval ION_BUILD_DIRECTORY=$(BUILD_DIRECTORY)$(INPUT_FILE_BASENAME)/)
	$(eval ION_OUTPUT_DIRECTORY=$(ION_BUILD_DIRECTORY)ions/)
	$(eval NEG_ION_BUILD_DIRECTORY=$(ION_BUILD_DIRECTORY)negative_build/)
	$(eval POS_ION_BUILD_DIRECTORY=$(ION_BUILD_DIRECTORY)positive_build/)
	$(eval ORCA_INPUT_BUILD_DIRECTORY=$(ION_BUILD_DIRECTORY)orca_inputs/)
	$(eval ORCA_OUTPUT_BUILD_DIRECTORY=$(ION_BUILD_DIRECTORY)orca_outputs/)
	$(eval MOBCAL_INPUT_BUILD_DIRECTORY=$(ION_BUILD_DIRECTORY)mobcal_inputs/)
	$(eval MOBCAL_OUTPUT_BUILD_DIRECTORY=$(ION_BUILD_DIRECTORY)mobcal_outputs/)
	$(eval PROCESS_IONS_OUTPUT=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/$(RESULTS_IONS_OUTPUT_DIRECTORY))
	$(eval PROCESS_CONFORMERS_OUTPUT=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/$(RESULTS_CONFORMERS_OUTPUT_DIRECTORY))
	$(eval PROCESS_ORCA_INPUTS_OUTPUT=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/$(RESULTS_ORCA_INPUTS_OUTPUT_DIRECTORY))
	$(eval PROCESS_ORCA_OUTPUTS_OUTPUT=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/$(RESULTS_ORCA_OUTPUTS_OUTPUT_DIRECTORY))
	$(eval PROCESS_MOBCAL_INPUTS_OUTPUT=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/$(RESULTS_MOBCAL_INPUTS_OUTPUT_DIRECTORY))
	$(eval PROCESS_MOBCAL_OUTPUTS_OUTPUT=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/$(RESULTS_MOBCAL_OUTPUTS_OUTPUT_DIRECTORY))
	$(eval INPUT_PATH=$(shell realpath "$(INPUT_FILE)"))
	$(eval TIMESTAMP_LOG=$(RESULTS_DIRECTORY)$(INPUT_FILE_BASENAME)/runtime.log)
	@printf "Cores=%s\nRAM per Core=%s\nStart Time=%s\n" $(AVAILABLE_CORES) $(AVAILABLE_RAM) `$(CURRENT_TIMESTAMP_COMMAND)` > $(TIMESTAMP_LOG)
	mkdir -p "$(ION_BUILD_DIRECTORY)"
	mkdir -p "$(ION_OUTPUT_DIRECTORY)"
	mkdir -p "$(NEG_ION_BUILD_DIRECTORY)"
	mkdir -p "$(POS_ION_BUILD_DIRECTORY)"
	mkdir -p "$(ORCA_OUTPUT_BUILD_DIRECTORY)"
	mkdir -p "$(MOBCAL_INPUT_BUILD_DIRECTORY)"
	mkdir -p "$(MOBCAL_OUTPUT_BUILD_DIRECTORY)"
	$(if $(wildcard $(PROCESS_IONS_OUTPUT)), rm -f "$(PROCESS_IONS_OUTPUT)"*, mkdir -p "$(PROCESS_IONS_OUTPUT)")
	$(if $(wildcard $(PROCESS_CONFORMERS_OUTPUT)), rm -f "$(PROCESS_CONFORMERS_OUTPUT)"*, mkdir -p "$(PROCESS_CONFORMERS_OUTPUT)")
	$(if $(wildcard $(PROCESS_ORCA_INPUTS_OUTPUT)), rm -f "$(PROCESS_ORCA_INPUTS_OUTPUT)"*, mkdir -p "$(PROCESS_ORCA_INPUTS_OUTPUT)")
	$(if $(wildcard $(PROCESS_ORCA_OUTPUTS_OUTPUT)), rm -f "$(PROCESS_ORCA_OUTPUTS_OUTPUT)"*, mkdir -p "$(PROCESS_ORCA_OUTPUTS_OUTPUT)")
	$(if $(wildcard $(PROCESS_MOBCAL_INPUTS_OUTPUT)), rm -f "$(PROCESS_MOBCAL_INPUTS_OUTPUT)"*, mkdir -p "$(PROCESS_MOBCAL_INPUTS_OUTPUT)")
	$(if $(wildcard $(PROCESS_MOBCAL_OUTPUTS_OUTPUT)), rm -f "$(PROCESS_MOBCAL_OUTPUTS_OUTPUT)"*, mkdir -p "$(PROCESS_MOBCAL_OUTPUTS_OUTPUT)")
	@printf "Setup Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	make --ignore-errors --makefile="$(CUR_MAKEFILE)" --directory="$(NEG_ION_BUILD_DIRECTORY)" generate-charged-particles INPUT_FILE="$(INPUT_PATH)" OUTPUT_DIRECTORY="$(PROCESS_IONS_OUTPUT)" CHARGE=-1
	@printf "Negative Particle Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	make --ignore-errors --makefile="$(CUR_MAKEFILE)" --directory="$(POS_ION_BUILD_DIRECTORY)" generate-charged-particles INPUT_FILE="$(INPUT_PATH)" OUTPUT_DIRECTORY="$(PROCESS_IONS_OUTPUT)" CHARGE=1
	@printf "Postive Particle Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	@make --no-print-directory remove-build-directory DIRECTORY=$(NEG_ION_BUILD_DIRECTORY)
	@make --no-print-directory remove-build-directory DIRECTORY=$(POS_ION_BUILD_DIRECTORY)
	@if [[ -z `ls $(PROCESS_IONS_OUTPUT)` ]]; then exit 1; fi
	@make --no-print-directory divider
	for FILE in `ls $(PROCESS_IONS_OUTPUT)*.xyz`; do MULTIPLICITY_STRING=`$(CHECK_MULTIPLICITY_STRING) $$FILE`; [[ -z $$MULTIPLICITY_STRING ]] && ION_MULTIPLICITY=1 || ION_MULTIPLICITY=2; CHECK_STRING=`$(CHECK_CHARGE_STRING) $$FILE`; [[ -z $$CHECK_STRING ]] && ION_CHARGE="1" || ION_CHARGE="-1"; CONFORMER_DIRECTORY="$(ION_OUTPUT_DIRECTORY)`basename $$FILE`/"; mkdir -p $$CONFORMER_DIRECTORY; make --no-print-directory divider; make --makefile="$(CUR_MAKEFILE)" --directory=$$CONFORMER_DIRECTORY generate-conformers SIMILARITY_THRESHOLD=$(PIPELINE_SIMILARITY_THRESHOLD) INPUT_FILE="$$FILE" OUTPUT_DIRECTORY="$(PROCESS_CONFORMERS_OUTPUT)" CHARGE=$$ION_CHARGE MULTIPLICITY=$$ION_MULTIPLICITY; done
	@printf "Conformer Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	for FILE in `find $(ION_OUTPUT_DIRECTORY) -name $(MASS_COMPARISON_OUTPUT_DIRECTORY)`; do CHECK_STRING=`$(CHECK_CHARGE_STRING) $$FILE`; MULTIPLICITY_STRING=`$(CHECK_MULTIPLICITY_STRING) $$FILE`; [[ -z $$CHECK_STRING ]] && ION_CHARGE="1" || ION_CHARGE="-1"; [[ -z $$MULTIPLICITY_STRING ]] && ION_MULTIPLICITY=1 || ION_MULTIPLICITY=2; make --no-print-directory divider; make --makefile="$(CUR_MAKEFILE)" --directory=$$FILE generate-orca-inputs CHARGE=$$ION_CHARGE MULTIPLICITY=$$ION_MULTIPLICITY INPUT_DIRECTORY="." OUTPUT_DIRECTORY="$(PROCESS_ORCA_INPUTS_OUTPUT)"; done
	@make --no-print-directory remove-build-directory DIRECTORY=$(ION_OUTPUT_DIRECTORY)
	@printf "ORCA Setup Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	make --makefile="$(CUR_MAKEFILE)" --directory="$(ORCA_OUTPUT_BUILD_DIRECTORY)" generate-orca-outputs INPUT_DIRECTORY="$(PROCESS_ORCA_INPUTS_OUTPUT)" OUTPUT_DIRECTORY="$(PROCESS_ORCA_OUTPUTS_OUTPUT)"
	@printf "ORCA Run Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory remove-build-directory DIRECTORY=$(ORCA_OUTPUT_BUILD_DIRECTORY)
	@make --no-print-directory divider
	make --makefile="$(CUR_MAKEFILE)" --directory="$(MOBCAL_INPUT_BUILD_DIRECTORY)" generate-mobcal-inputs INPUT_DIRECTORY="$(PROCESS_ORCA_OUTPUTS_OUTPUT)" OUTPUT_DIRECTORY="$(PROCESS_MOBCAL_INPUTS_OUTPUT)" CHARGE_SETTING=$(PIPELINE_MFJ_CHARGE_SETTING) MFJ_CYCLES=$(PIPELINE_MFJ_CYCLES) MFJ_VELOCITY_INTEGRATION=$(PIPELINE_MFJ_VELOCITY_INTEGRATION) MFJ_IMPACT_INTEGRATION=$(PIPELINE_MFJ_IMPACT_INTEGRATION) MFJ_CARRIER_GAS=$(PIPELINE_MFJ_CARRIER_GAS) MFJ_TEMPERATURE=$(PIPELINE_MFJ_TEMPERATURE) MFJ_EMPIRICAL_CORRECTION=`$(EMPIRICAL_CORRECTIONS_CHECK) $(PROCESS_ORCA_OUTPUTS_OUTPUT)*.out`
	@make --no-print-directory remove-build-directory DIRECTORY=$(MOBCAL_INPUT_BUILD_DIRECTORY)
	@printf "MobCal-MPI Setup Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	make --makefile="$(CUR_MAKEFILE)" --directory="$(MOBCAL_OUTPUT_BUILD_DIRECTORY)" generate-mobcal-outputs INPUT_DIRECTORY="$(PROCESS_MOBCAL_INPUTS_OUTPUT)" OUTPUT_DIRECTORY="$(PROCESS_MOBCAL_OUTPUTS_OUTPUT)"
	@make --no-print-directory remove-build-directory DIRECTORY=$(MOBCAL_OUTPUT_BUILD_DIRECTORY)
	@printf "MobCal-MPI Run Complete Time=%s\n" `$(CURRENT_TIMESTAMP_COMMAND)` >> $(TIMESTAMP_LOG)
	@make --no-print-directory divider
	@printf "Processing complete.\nResults stored in '%s'\n" "$(RESULTS_DIRECTORY)"

generate-charged-particles:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE)" INPUT_VARIABLE="CHARGE"
	$(eval CREST_OPERATION=$(shell make --no-print-directory --makefile="$(CUR_MAKEFILE)" charge-to-conjugate CHARGE="$(CHARGE)"))
	$(eval LOG_FILE=$(if $(OUTPUT_DIRECTORY),$(shell pwd),$(OUTPUT_DIRECTORY))/$(if $(INPUT_FILE),$(shell printf "%s.crest%slog" $(shell basename "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*//") "$(CREST_OPERATION)"),tester.log))
	$(eval NEUTRAL_FILE_NAME=$(shell pwd)/$(notdir $(INPUT_FILE:.sdf=.xyz)))
	$(if $(shell $(XYZ_FILE_FORMAT_CHECKER) "$(INPUT_FILE)" 2> "$(LOG_FILE)" | grep "XYZ"), cp "$(INPUT_FILE)" "$(NEUTRAL_FILE_NAME)", $(OPEN_BABEL) -isdf "$(INPUT_FILE)" -oxyz --gen3D > "$(NEUTRAL_FILE_NAME)")
	$(eval ION_BUILD=$(shell pwd)/ion_build)
	mkdir -p "$(ION_BUILD)"
	make --ignore-errors --makefile="$(CUR_MAKEFILE)" --directory=$(ION_BUILD) generate-ion INPUT_FILE="$(NEUTRAL_FILE_NAME)" OUTPUT_DIRECTORY="$(realpath $(OUTPUT_DIRECTORY))" CHARGE="$(CHARGE)"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" divider
	$(eval CONJUGATE_BUILD=$(shell pwd)/conjugate_build)
	mkdir -p "$(CONJUGATE_BUILD)"
	make --makefile="$(CUR_MAKEFILE)" --directory=$(CONJUGATE_BUILD) generate-conjugate INPUT_FILE="`realpath \"$(INPUT_FILE)\"`" OUTPUT_DIRECTORY="$(realpath $(OUTPUT_DIRECTORY))" CHARGE="$(CHARGE)"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" divider
	$(eval DIMER_BUILD=$(shell pwd)/conjugate_dimer_build)
	mkdir -p "$(DIMER_BUILD)"
	$(eval OUTPUT_FILE=$(OUTPUT_DIRECTORY)/$(notdir $(shell printf "%s" "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*/_$(CREST_OPERATION)\.xyz/")))
	make --makefile="$(CUR_MAKEFILE)" --directory=$(DIMER_BUILD) generate-dimer INPUT_FILE="`realpath \"$(OUTPUT_FILE)\"`" NEUTRAL_FILE="$(NEUTRAL_FILE_NAME)" OUTPUT_DIRECTORY="$(realpath $(OUTPUT_DIRECTORY))" CHARGE=$(if $(shell printf "%s" $(CREST_OPERATION) | grep "deprotonate"),-1,1)

generate-conjugate:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE)" INPUT_VARIABLE="CHARGE"
	$(eval CREST_OPERATION=$(shell make --no-print-directory --makefile="$(CUR_MAKEFILE)" charge-to-conjugate CHARGE="$(CHARGE)"))
	$(eval LOG_FILE=$(if $(OUTPUT_DIRECTORY),$(shell pwd),$(OUTPUT_DIRECTORY))/$(if $(INPUT_FILE),$(shell printf "%s.crest%slog" $(shell basename "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*//") "$(CREST_OPERATION)"),tester.log))
	@printf "Determining ions of '%s' with %s operation of crest utility\n" "$(INPUT_FILE)" "$(CREST_OPERATION)"
	crest "$(INPUT_FILE)" --$(CREST_OPERATION) --T $(AVAILABLE_CORES) &> "$(LOG_FILE)"
	@if ! [[ -f $(CREST_OPERATION)d.xyz ]]; then printf "Unable to find conformers for '%s' with %s operation.\n" "$(INPUT_FILE)" "$(CREST_OPERATION)"; exit 1; fi
	crest --splitfile $(CREST_OPERATION)d.xyz 1 1 &>> "$(LOG_FILE)"
	$(eval OUTPUT_FILE=$(OUTPUT_DIRECTORY)/$(notdir $(shell printf "%s" "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*/_$(CREST_OPERATION)\.xyz/")))
	cp "`find SPLIT/ -name struc.xyz`" "$(OUTPUT_FILE)"

generate-ion:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE)" INPUT_VARIABLE="CHARGE"
	$(eval ION_TYPE=$(shell make --no-print-directory --makefile="$(CUR_MAKEFILE)" charge-to-ion CHARGE="$(CHARGE)"))
	$(eval LOG_FILE=$(if $(OUTPUT_DIRECTORY),$(OUTPUT_DIRECTORY),$(shell pwd))/$(if $(INPUT_FILE),$(shell printf "%s.crest%slog" $(shell basename "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*//") $(ION_TYPE)),tester.log))
	@printf "Determining %s of '%s' with crest utility\n" $(ION_TYPE) "$(INPUT_FILE)"
	crest "$(INPUT_FILE)" --chrg $(CHARGE) --uhf 1 --gfnff --T $(AVAILABLE_CORES) &> "$(LOG_FILE)"
	@if ! [[ -f crest_best.xyz ]]; then printf "Unable to find %s for '%s' with crest.\n" "$(ION_TYPE)" "$(INPUT_FILE)"; exit 1; fi
	cp crest_best.xyz "$(OUTPUT_DIRECTORY)/$(notdir $(INPUT_FILE:.xyz=_$(ION_TYPE).xyz))"

generate-dimer:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(NEUTRAL_FILE)" INPUT_VARIABLE="NEUTRAL_FILE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE)" INPUT_VARIABLE="CHARGE"
	$(eval DIMER_FILE=$(if $(shell printf "%s" $(CHARGE) | grep "-"),neg, pos)_dimer.xyz)
	$(eval LOG_FILE=$(OUTPUT_DIRECTORY)/$(shell printf "%s_%s.dimerlog" $(shell basename "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*//") $(if $(shell printf "%s" $(CHARGE) | grep "-"),neg, pos)))
	$(DIMER_APPEND) "$(INPUT_FILE)" "$(NEUTRAL_FILE)" "$(DIMER_FILE)"
	crest "$(DIMER_FILE)" --chrg $(CHARGE) --gfnff --wscal 1.0 --nci --noreftopo --T $(AVAILABLE_CORES) &> "$(LOG_FILE)"
	@if ! [[ -f crest_best.xyz ]]; then printf "Unable to find dimer for '%s' with crest.\n" "$(INPUT_FILE)"; exit 1; fi
	cp crest_best.xyz "$(OUTPUT_DIRECTORY)/$(notdir $(INPUT_FILE:.xyz=_dimer.xyz))"

generate-conformers:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE)" INPUT_VARIABLE="CHARGE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(SIMILARITY_THRESHOLD)" INPUT_VARIABLE="SIMILARITY_THRESHOLD"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MULTIPLICITY)" INPUT_VARIABLE="MULTIPLICITY"
	$(eval COSINE_SIMILARITY=$(if $(shell $(CHECK_COSINE_SIMILAIRITY) $(MULTIPLICITY)),--uhf 1 ,)) 
	$(eval LOG_FILE=$(OUTPUT_DIRECTORY)/$(shell printf "%s.conformerlog" $(shell basename "$(INPUT_FILE)" | sed "s/ /_/g" | sed "s/\..*//")))
	@printf "Generating conformers with crest utility for '%s'\n" "$(INPUT_FILE)"
	crest "$(INPUT_FILE)" --gfn2 --chrg $(CHARGE) --T $(AVAILABLE_CORES) $(COSINE_SIMILARITY)&> "$(LOG_FILE)"
	@printf "Determining unique conformers based upon a %s%% threshold\n" $(SIMILARITY_THRESHOLD)
	$(ION_SPLITTER) . $(CONFORMERS_FILE) $(basename $(notdir $(INPUT_FILE))) $(ION_SPLITTER_OUTPUT_DIRECTORY)
	$(COMPARE_MASSES) $(ION_SPLITTER_OUTPUT_DIRECTORY) $(SIMILARITY_THRESHOLD)
	cp $(ION_SPLITTER_OUTPUT_DIRECTORY)/$(MASS_COMPARISON_OUTPUT_DIRECTORY)/*.gjf "$(OUTPUT_DIRECTORY)"

generate-orca-inputs:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_DIRECTORY)" INPUT_VARIABLE="INPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE)" INPUT_VARIABLE="CHARGE"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MULTIPLICITY)" INPUT_VARIABLE="MULTIPLICITY"
	$(eval ORCA_INPUTS_DIRECTORY=ORCAInputs)
	cp -n $(INPUT_DIRECTORY)/*.gjf .
	$(GJF_ORCA_CONVERSION) . $(AVAILABLE_RAM) $(AVAILABLE_CORES) $(CHARGE) $(MULTIPLICITY)
	cp -n $(ORCA_INPUTS_DIRECTORY)/* "$(OUTPUT_DIRECTORY)"

generate-orca-outputs:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_DIRECTORY)" INPUT_VARIABLE="INPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@printf "Processing all files in '%s' with orca utility.\n" "$(INPUT_DIRECTORY)"
	@printf "Copying all files from input directory.\n"
	$(shell cp -n $(INPUT_DIRECTORY)/*.inp .)
	$(foreach file,$(notdir $(wildcard *.inp)),printf "Processing '$(file)' with orca utility.\n";$(ORCA) "$(file)" &> "$(file:.inp=.out)";)
	rm -f *.scfhess.out
	cp -n *.out "$(OUTPUT_DIRECTORY)"
	@printf "Extracting thermochem information.\n"
	$(THERMO_CHEM_EXTRACTOR) . "$(OUTPUT_DIRECTORY)"

generate-mobcal-inputs:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_DIRECTORY)" INPUT_VARIABLE="INPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHARGE_SETTING)" INPUT_VARIABLE="CHARGE_SETTING"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MFJ_CYCLES)" INPUT_VARIABLE="MFJ_CYCLES"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MFJ_VELOCITY_INTEGRATION)" INPUT_VARIABLE="MFJ_VELOCITY_INTEGRATION"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MFJ_IMPACT_INTEGRATION)" INPUT_VARIABLE="MFJ_IMPACT_INTEGRATION"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MFJ_CARRIER_GAS)" INPUT_VARIABLE="MFJ_CARRIER_GAS"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MFJ_EMPIRICAL_CORRECTION)" INPUT_VARIABLE="MFJ_EMPIRICAL_CORRECTION"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(MFJ_TEMPERATURE)" INPUT_VARIABLE="MFJ_TEMPERATURE"
	$(eval CREATOR_LOG=$(realpath $(OUTPUT_DIRECTORY))/$(MFJ_CREATOR_LOG))
	@printf "Preparing MobCal-MPI inputs.\n"
	cp -n $(INPUT_DIRECTORY)/*.out .
	$(eval MFJ_CREATOR_STRING=$(MFJ_CYCLES),$(MFJ_VELOCITY_INTEGRATION),$(MFJ_IMPACT_INTEGRATION),$(MFJ_CARRIER_GAS),$(MFJ_EMPIRICAL_CORRECTION),$(MFJ_TEMPERATURE))
	$(MFJ_CREATOR) "." "$(SDF2TINKERXYZ)" "$(CHARGE_SETTING)" "$(MFJ_CREATOR_STRING)" >> $(CREATOR_LOG)
	@SUCCESSRATE=(`cat $(CREATOR_LOG) | grep "Success Rate:"`); if ! [[ $${SUCCESSRATE[-1]} -gt 0 ]]; then printf "Unable to process any ORCA output files.\n"; exit 1; fi
	cp $(MFJ_CREATOR_OUTPUTS_DIRECTORY)/*.mfj "$(OUTPUT_DIRECTORY)"

generate-mobcal-outputs: $(SCRIPTS_DIRECTORY)mobcal_mpi.exe
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_DIRECTORY)" INPUT_VARIABLE="INPUT_DIRECTORY"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@printf "Processing all files in '%s' with mpirun utility.\n" "$(INPUT_DIRECTORY)"
	$(shell cp -n $(INPUT_DIRECTORY)/*.mfj .)
	$(foreach file,$(wildcard *.mfj),printf "Processing '$(file)' with mpirun utility.\n"; $(MPIRUN) -np $(AVAILABLE_CORES) "$<" "$(file)" "$(file:.mfj=.mout)";)
	cp -n *.mout $(OUTPUT_DIRECTORY)

$(SCRIPTS_DIRECTORY)mobcal_mpi.exe: $(wildcard $(FORTRAN_DIRECTORY)*.f)
	$(FORTRAN_COMPILER) -o $@ $< $(FORTRAN_FLAGS)

convert-smiles: update
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHEMICAL_STRING)" INPUT_VARIABLE="CHEMICAL_STRING"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHEMICAL_NAME)" INPUT_VARIABLE="CHEMICAL_NAME"
	make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-chemical-exists CHEMICAL_NAME="$(CHEMICAL_NAME)"
	@printf "%s\n" "$(CHEMICAL_STRING)" > $(CHEMICAL_NAME).smi; $(OPEN_BABEL) -i SMI "$(CHEMICAL_NAME).smi" -o SDF -O "$(CHEMICAL_NAME).sdf"
	@if [[ -z `cat "$(CHEMICAL_NAME).sdf"` ]]; then printf "Open Babel processing error. Check inputs.\n"; exit 1; fi
	make --makefile="$(CUR_MAKEFILE)" add-chemical INPUT_FILE="$(CHEMICAL_NAME).sdf" CHEMICAL_NAME="$(CHEMICAL_NAME)"

download-chemical: update
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHEMICAL_CID)" INPUT_VARIABLE="CHEMICAL_CID"
	make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-chemical-exists CHEMICAL_NAME="$(CHEMICAL_CID)"
	$(eval CHEMICAL_FILE_NAME=COMPOUND_CID_$(CHEMICAL_CID))
	$(eval CHEMICAL_FILE_PATH=COMPOUND_CID_$(CHEMICAL_CID).sdf)
	$(eval LOC_PUBCHEM_DATA=data.csv)
	@sleep 1s # Sleep by one second to adhere to PubChem PUG REST usage policy
	@wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(CHEMICAL_CID)/record/SDF?record_type=3d -O "$(CHEMICAL_FILE_PATH)" &>> $(MAKEFILE_LOG)
	@wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(CHEMICAL_CID)/property/IUPACName\,MolecularFormula\,ConnectivitySMILES\,MolecularWeight\,Volume3D/CSV -O $(LOC_PUBCHEM_DATA) &>> $(MAKEFILE_LOG)
	make --no-print-directory --makefile="$(CUR_MAKEFILE)" add-chemical CHEMICAL_NAME=$(CHEMICAL_FILE_NAME) INPUT_FILE=$(CHEMICAL_FILE_PATH)
	make --no-print-directory --makefile="$(CUR_MAKEFILE)" process-pubchem-data PUBCHEM_DATA="$(LOC_PUBCHEM_DATA)" OUTPUT_DIRECTORY="."
	make --no-print-directory --makefile="$(CUR_MAKEFILE)" add-chemical CHEMICAL_NAME=$(CHEMICAL_FILE_NAME) INPUT_FILE=$(PUBCHEM_DATA_OUTPUT_FILE)
	@printf "Downloaded CID '%s' to '%s'\n" "$(CHEMICAL_CID)" "$(CHEMICAL_FILE_NAME)"

add-chemical: update
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHEMICAL_NAME)" INPUT_VARIABLE="CHEMICAL_NAME" 
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(INPUT_FILE)" INPUT_VARIABLE="INPUT_FILE"
	$(eval CHEMICAL_NAME_DIRECTORY=$(RESULTS_DIRECTORY)$(CHEMICAL_NAME)/$(RESULTS_CHEMICAL_STRUCTURE_DIRECTORY))
	mkdir -p "$(CHEMICAL_NAME_DIRECTORY)"
	mv "$(INPUT_FILE)" "$(CHEMICAL_NAME_DIRECTORY)/"
	@printf "Moved '%s' to '%s'\n" "$(INPUT_FILE)" "$(CHEMICAL_NAME_DIRECTORY)"

process-pubchem-data:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(PUBCHEM_DATA)" INPUT_VARIABLE="PUBCHEM_DATA"
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(OUTPUT_DIRECTORY)" INPUT_VARIABLE="OUTPUT_DIRECTORY"
	@if ! [[ -f "$(PUBCHEM_DATA)" ]]; then printf "Unable to find file '%s'\n" "$(PUBCHEM_DATA)"; exit 1; fi
	$(eval OUTPUT_DATA_PATH=$(OUTPUT_DIRECTORY)/$(PUBCHEM_DATA_OUTPUT_FILE))
	@printf "" > "$(OUTPUT_DATA_PATH)"
	@ION_KEYS=(`head -n 1 $(PUBCHEM_DATA) | sed "s/,/ /g"`); ION_VALS=(`tail -n 1 "$(PUBCHEM_DATA)" | sed "s/ /_/g" | sed "s/,/ /g"`); for INDEX in `seq 1 $$(( $${#ION_KEYS[@]} - 1 ))`; do printf "%s: %s\n" $${ION_KEYS[INDEX]} "$${ION_VALS[INDEX]}" | sed "s/\"//g" | sed "s/_/ /g" >> "$(OUTPUT_DATA_PATH)"; done

check-chemical-exists:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(CHEMICAL_NAME)" INPUT_VARIABLE="CHEMICAL_NAME" 
	@if [[ -d "$(RESULTS_CHEMICAL_STRUCTURE_DIRECTORY)$(CHEMICAL_NAME)" ]]; then printf "ERROR: Chemical Directory for '%s' exists. Unable to add chemical.\n" "$(CHEMICAL_NAME)"; exit 2; fi

check-input: 
	@if [[ -z "$(INPUT_VALUE)" ]]; then printf "ERROR: Invalid value provided for '%s' input.\n" "$(INPUT_VARIABLE)"; exit 1; fi

charge-to-conjugate:
	$(eval CREST_OPERATION=$(if $(findstring -1,$(CHARGE)),deprotonate,protonate))
	@printf "%s\n" $(CREST_OPERATION)

charge-to-ion:
	$(eval CREST_OPERATION=$(if $(findstring -1,$(CHARGE)),radical-anion,radical-cation))
	@printf "%s\n" $(CREST_OPERATION)

remove-build-directory:
	@make --no-print-directory --makefile="$(CUR_MAKEFILE)" check-input INPUT_VALUE="$(DIRECTORY)" INPUT_VARIABLE="DIRECTORY"
	$(if $(CLEANUP_BUILD),rm -rf $(DIRECTORY),exit 0)

divider:
	@printf "\n"; for ITEM in `seq 1 $(DIVIDER_COUNT)`; do printf "%s" "-"; done; printf "\n\n"

.PHONY: clean
clean:
	rm -r "$(BUILD_DIRECTORY)"

