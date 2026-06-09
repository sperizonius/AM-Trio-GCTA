## Script to run the different steps for Design 1-4 in one go, it is also possible to parallelize these jobs if your server allows this
## Sofie Perizonius - 9-6-2026

#!/bin/bash

set -e

echo "[$(date)] Starting full pipeline for DESIGN 1–4..."

BASE_DIR="/data/sofie/Internship_2/simulations"

#for i in 1
for i in 1 2 3 4
do
  echo "----------------------------------------"
  echo "[$(date)] Starting DESIGN ${i}"
  echo "----------------------------------------"

  # Step 1: Go to simulation data folder
  cd ${BASE_DIR}/simulated_data/DESIGN_${i}

  echo "[$(date)] Running Python simulation for DESIGN ${i}..."
  python3 ../../simulation_scripts/design_${i}.py > ../../output_error_files/design_${i}_sim.o 2> ../../output_error_files/design_${i}_sim.e 

  echo "[$(date)] Running PLINK step for DESIGN ${i}..."
  bash ../../simulation_scripts/plink_design_${i}.sh > ../../output_error_files/design_${i}_plink.o 2> ../../output_error_files/design_${i}_plink.e 

  # Step 2: Go to main simulations folder
  cd ${BASE_DIR}

  echo "[$(date)] Running Julia analysis for DESIGN ${i}..."
  julia simulation_scripts/trio_design_${i}_standard.jl > output_error_files/design_${i}_standard_analysis.o 2> output_error_files/design_${i}_standard_analysis.e &
  
  PID1=$!
  echo "PID standard mode design ${i} ${PID1}"
  
  julia simulation_scripts/trio_design_${i}_simple_adapted.jl > output_error_files/design_${i}_simple_adapted_analysis.o 2> output_error_files/design_${i}_simple_adapted_analysis.e &
  
  PID2=$!
  echo "PID adapted mode design ${i} ${PID2}"
  
  # wait for both jobs to finish
  wait $PID1
  wait $PID2

  echo "[$(date)] Finished DESIGN ${i}"

done

echo "========================================"
echo "[$(date)] DESIGN 1-4 FINISHED!"
echo "========================================"
