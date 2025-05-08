#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=5:mem=20gb
#PBS -J 0-2

# Load required modules
module load singularity

# Change to the directory the script was submitted from (this is inside root/build/Release)
cd $PBS_O_WORKDIR

#echo the current working directory
echo "Current working directory: $(pwd)"

# Define parameter arrays
KAPPA_VALUES=(0 0.001 0.01 0.02 0.03 0.04 0.05)
STRAIN_TYPES=("diastolic" "systolic" "no_strain")

# Calculate which parameters to use based on PBS array index
# Each strain type has 7 kappa values
NUM_KAPPAS=${#KAPPA_VALUES[@]}
STRAIN_INDEX=$(( PBS_ARRAY_INDEX / NUM_KAPPAS ))
KAPPA_INDEX=$(( PBS_ARRAY_INDEX % NUM_KAPPAS ))

# Get the specific parameter values for this job
STRAIN_TYPE=${STRAIN_TYPES[$STRAIN_INDEX]}
KAPPA=${KAPPA_VALUES[$KAPPA_INDEX]}

# Use PBS_ARRAY_INDEX as the seed
SEED=${PBS_ARRAY_INDEX}

# Set other parameters with default values
ANGLE=0.01
SHIFT_BLOCK="false"
D_ECS=2.5
D_ICS=1.0
CORES=5
IS_DEFORMED="true"
STRAIN_STEP_SIZE=100
PARTICLES=10

# For no_strain case, set isDeformed to false for better performance
if [ "$STRAIN_TYPE" == "no_strain" ]; then
    IS_DEFORMED="false"
fi

# Pull Singularity image if file is not there (you have to do it before running the script)
# Command: singularity pull simulation.sif docker://ignasiialemany/mcrw-cdti:latest
if [ ! -f "simulation.sif" ]; then
    singularity pull simulation.sif docker://ignasiialemany/mcrw-cdti:latest
fi

# Run simulation using Singularity
singularity exec \
    ./simulation.sif \
    /app/build/Release/3DRandomWalk \
        --sequence-file sequence.yaml \
        --geometry-file geometry_1.mat \
        --strain-type ${STRAIN_TYPE} \
        --particles ${PARTICLES} \
        --kappa ${KAPPA} \
        --seed ${SEED} \
        --angle ${ANGLE} \
        --shift-block ${SHIFT_BLOCK} \
        --Decs ${D_ECS} \
        --Dics ${D_ICS} \
        --cores ${CORES} \
        --deformed ${IS_DEFORMED} \
        --strain-step-size ${STRAIN_STEP_SIZE} \
        --job-id ${PBS_ARRAY_INDEX}
