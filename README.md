## Reduced_Order_Slosh_Estimation_Model

# System Dependencies (Fedora)
OpenFOAM requires a full C++ + MPI + scientific toolchain.

Install everything with:

```bash
sudo dnf groupinstall "Development Tools"

sudo dnf install \
    gcc gcc-c++ make cmake \
    git \
    flex bison \
    zlib-devel \
    boost-devel \
    openmpi openmpi-devel \
    scotch scotch-devel \
    paraview \
    qt5-qtbase-devel qt5-qtsvg-devel qt5-qttools-devel \
    libXt-devel \
    gperftools-devel \
    fftw-devel
