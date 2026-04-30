#!/bin/bash

echo "Loading OpenFOAM 11 environment..."

export PROJECT_HOME=$HOME/Reduced_Order_Slosh_Estimation_Model

export FOAM_INST_DIR=$PROJECT_HOME

# Use system OpenMPI (same as before)
export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

# Source OpenFOAM 11
source $PROJECT_HOME/OpenFOAM-11/etc/bashrc

echo "OpenFOAM version:"
foamVersion

echo "which interFoam:"
which interFoam
