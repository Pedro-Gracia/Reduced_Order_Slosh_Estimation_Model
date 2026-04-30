#!/bin/bash

# ============================================================
# OpenFOAM Case Cleaner
# Resets a case to pre-simulation state
# Usage:
#   clean_case.sh <case_path>
# Example:
#   clean_case.sh sims/SIM_Cavity
# ============================================================

# ----------- Input handling -----------
CASE_DIR="$1"

# Remove processor directories
rm -rf processor*

if [ -z "$CASE_DIR" ]; then
    echo "ERROR: No case path provided"
    echo "Usage: clean_case.sh <case_path>"
    return 1 2>/dev/null || exit 1
fi

if [ ! -d "$CASE_DIR" ]; then
    echo "ERROR: Case directory does not exist: $CASE_DIR"
    return 1 2>/dev/null || exit 1
fi

echo "----------------------------------------"
echo "Cleaning OpenFOAM case:"
echo "$CASE_DIR"
echo "----------------------------------------"

cd "$CASE_DIR" || exit 1

# ----------- Remove time directories -----------
echo "Removing time directories..."

# Keep only "0", remove everything else that looks like a time folder
for d in $(ls -d [0-9]* 2>/dev/null); do
    if [ "$d" != "0" ]; then
        echo "Deleting time directory: $d"
        rm -rf "$d"
    fi
done

# ----------- Remove postProcessing -----------
if [ -d "postProcessing" ]; then
    echo "Removing postProcessing/"
    rm -rf postProcessing
fi

# ----------- Remove logs -----------
echo "Removing log files..."
rm -f log.*
rm -f *.log

# ----------- Remove processor dirs (parallel runs) -----------
echo "Removing processor directories..."
rm -rf processor*

# ----------- Optional: remove mesh -----------
read -p "Remove mesh (constant/polyMesh)? [y/N]: " REMOVE_MESH

if [[ "$REMOVE_MESH" == "y" || "$REMOVE_MESH" == "Y" ]]; then
    echo "Removing mesh..."
    rm -rf constant/polyMesh
else
    echo "Keeping mesh."
fi

# ----------------------------------------
# Reset alpha.water 
# ----------------------------------------
echo "Resetting 0/alpha.water..."

cat > 0/alpha.water << EOF
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.water;
}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    walls
    {
        type zeroGradient;
    }
}
EOF

# --------------------------------------------------
# Remove logs
# --------------------------------------------------
rm -f log.*

echo "----------------------------------------"
echo "Case cleaned successfully"
echo "----------------------------------------"
