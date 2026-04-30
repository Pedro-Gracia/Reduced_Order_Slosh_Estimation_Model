#!/bin/bash

set -e

# ============================================
# Check input
# ============================================
if [ -z "$1" ]; then
    echo "Usage: $0 <relative_case_path>"
    exit 1
fi

CASE_DIR=$1

if [ ! -d "$CASE_DIR" ]; then
    echo "Error: Directory '$CASE_DIR' does not exist"
    exit 1
fi

echo "-----------------------------------------"
echo "Setting up liquidCOM in:"
echo "$CASE_DIR"
echo "-----------------------------------------"

cd "$CASE_DIR"

DIR="liquidCOMFunction"
mkdir -p $DIR/Make
cd $DIR

# ============================================
# liquidCOM.H
# ============================================
cat > liquidCOM.H << 'EOF'
#ifndef liquidCOM_H
#define liquidCOM_H

#include "functionObject.H"
#include "volFields.H"
#include <fstream>

namespace Foam
{
namespace functionObjects
{

class liquidCOM
:
    public functionObject
{
    const fvMesh& mesh_;
    std::ofstream file_;

    scalar xCOM0_;
    bool initialized_;

public:

    TypeName("liquidCOM");

    liquidCOM
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual bool execute();
    virtual bool write();
    virtual wordList fields() const;
};

}
}

#endif
EOF

# ============================================
# liquidCOM.C
# ============================================
cat > liquidCOM.C << 'EOF'
#include "liquidCOM.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "polyMesh.H"

namespace Foam
{
namespace functionObjects
{

defineTypeNameAndDebug(liquidCOM, 0);
addToRunTimeSelectionTable(functionObject, liquidCOM, dictionary);

// ============================================================
// Constructor
// ============================================================

liquidCOM::liquidCOM
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name, runTime),
    mesh_
    (
        refCast<const fvMesh>
        (
            runTime.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault<word>("region", "region0")
            )
        )
    ),
    xCOM0_(0.0),
    initialized_(false)
{
    if (Pstream::master())
    {
        fileName outDir = runTime.globalPath()
                        / "postProcessing"
                        / "liquidCOM"
                        / "0";

        mkDir(outDir);

        fileName outFile = outDir / "liquidCOM.dat";

        if (isFile(outFile))
        {
            file_.open(outFile.c_str(), std::ios::app);
        }
        else
        {
            file_.open(outFile.c_str());
            file_ << "# Time xCOM zCOM volume q" << std::endl;
        }
    }
}

// ============================================================
// Execute (MAIN LOGIC)
// ============================================================

bool liquidCOM::execute()
{
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>("alpha.water");

    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();

    scalar sumAlphaVol  = 0.0;
    scalar sumAlphaXVol = 0.0;
    scalar sumAlphaZVol = 0.0;

    forAll(alpha, cellI)
    {
        const scalar a   = alpha[cellI];
        const scalar vol = V[cellI];
        const scalar x   = C[cellI].x();
        const scalar z   = C[cellI].z();

        const scalar m = a * vol;

        sumAlphaVol  += m;
        sumAlphaXVol += m * x;
        sumAlphaZVol += m * z;
    }

    reduce(sumAlphaVol,  sumOp<scalar>());
    reduce(sumAlphaXVol, sumOp<scalar>());
    reduce(sumAlphaZVol, sumOp<scalar>());

    scalar xCOM = 0.0;
    scalar zCOM = 0.0;

    if (sumAlphaVol > SMALL)
    {
        xCOM = sumAlphaXVol / sumAlphaVol;
        zCOM = sumAlphaZVol / sumAlphaVol;

        if (!initialized_)
        {
            xCOM0_ = xCOM;
            initialized_ = true;
        }
    }

    scalar q = xCOM - xCOM0_;

    if (Pstream::master() && file_.is_open())
    {
        file_ << mesh_.time().value() << " "
              << xCOM << " "
              << zCOM << " "
              << sumAlphaVol << " "
              << q << std::endl;
    }

    return true;
}

// ============================================================
// Write (unused)
// ============================================================

bool liquidCOM::write()
{
    return true;
}

// ============================================================
// Fields
// ============================================================

wordList liquidCOM::fields() const
{
    return wordList(1, "alpha.water");
}

}
}
EOF

# ============================================
# Make/files
# ============================================
cat > Make/files << 'EOF'
liquidCOM.C

LIB = $(FOAM_USER_LIBBIN)/libliquidCOM
EOF

# ============================================
# Make/options
# ============================================
cat > Make/options << 'EOF'
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude

LIB_LIBS = \
    -lfiniteVolume
EOF

# ============================================
# Compile
# ============================================
echo "Compiling liquidCOM..."
wclean
wmake libso

echo "-----------------------------------------"
echo "DONE"
echo "-----------------------------------------"