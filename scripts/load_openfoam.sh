#!/bin/bash

REPO_DIR="$HOME/Reduced_Order_Slosh_Estimation_Model"

# ------------------------------
# 1) System tools (CRITICAL)
# ------------------------------
export PATH=/usr/bin:/bin

# ------------------------------
# 2) MPI (must come BEFORE OpenFOAM)
# ------------------------------
export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

# ------------------------------
# 3) SCOTCH (your manual install)
# ------------------------------
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# ------------------------------
# 4) Load OpenFOAM
# ------------------------------
cd "$REPO_DIR/OpenFOAM-10"
source etc/bashrc

# ------------------------------
# 5) (Optional sanity check)
# ------------------------------
echo "OpenFOAM loaded:"
which wmake
which mpicc
