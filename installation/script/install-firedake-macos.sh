#!/bin/bash

set -e

PETSC_OPTIONS_BASE="
  --with-c2html=0 \
  --with-debugging=0 \
  --with-fortran-bindings=0 \
  --with-shared-libraries=1 \
  --with-strict-petscerrorcode \
  --download-cmake \
  --download-bison \
  --download-fftw \
  --download-mumps-avoid-mpi-in-place \
  --with-hdf5-dir=/opt/homebrew \
  --with-hwloc-dir=/opt/homebrew \
  --download-metis \
  --download-mumps \
  --download-netcdf \
  --download-pnetcdf \
  --download-ptscotch \
  --download-scalapack \
  --download-suitesparse \
  --download-superlu_dist \
  --download-slepc \
  --with-zlib \
"

#  --download-pastix \

 
PETSC_EXTRA_COMMON=" \
    --download-hpddm --download-libpng \
    --download-ctetgen --download-tetgen --download-triangle \
    --download-mmg --download-parmmg --download-p4est"

PETSC_EXTRA_REAL="--download-eigen --download-hypre --download-pragmatic"
PETSC_EXTRA_COMPLEX="--with-scalar-type=complex"
PETSC_INT64="--with-64-bit-indices"

function get_petsc_arch() {
    ARCH=$1
    echo arch-firedrake-${ARCH}
}

function get_env_name() {
    ARCH=$1
    echo venv-firedrake-${ARCH}
}

function build_petsc() {
    ARCH=$1
    PETSC_ARCH=$(get_petsc_arch ${ARCH})
    shift
    echo PETSC_ARCH=$PETSC_ARCH
    pushd $PETSC_DIR
    ./configure PETSC_ARCH=$PETSC_ARCH \
               --COPTFLAGS='-O3 -march=native -mtune=native' \
               --CXXOPTFLAGS='-O3 -march=native -mtune=native' \
               --FOPTFLAGS='-O3 -mtune=native' \
               $PETSC_OPTIONS_BASE $PETSC_EXTRA_COMMON "$@"
    make PETSC_ARCH=$PETSC_ARCH all
    make PETSC_ARCH=$PETSC_ARCH check
    popd
}

function install_firedrake() {
    ARCH=$1
    ENV_NAME=$(get_env_name ${ARCH})
    PETSC_ARCH=$(get_petsc_arch ${ARCH})
    SLEPC_DIR=$PETSC_DIR/$PETSC_ARCH

    python3 -m pip cache purge
    python3 -m venv --copies $ENV_NAME
    source $ENV_NAME/bin/activate
    mkdir -p $ENV_NAME/src

    echo ENV: CC=mpicc CXX=mpicxx
    echo ENV: PETSC_DIR=$PETSC_DIR
    echo ENV: PETSC_ARCH=$PETSC_ARCH
    echo ENV: SLEPC_DIR=$SLEPC_DIR 
    echo ENV: HDF5_MPI=ON HDF5_DIR=/opt/homebrew
    CC=mpicc CXX=mpicxx \
        PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH SLEPC_DIR=$SLEPC_DIR \
        HDF5_MPI=ON HDF5_DIR=/opt/homebrew \
        pip install -v --no-binary h5py --src $ENV_NAME/src \
            --editable "git+https://github.com/firedrakeproject/firedrake.git@release#egg=firedrake[check,vtk,netgen,slepc]"
        # pip install --no-binary h5py "firedrake @ git+https://github.com/firedrakeproject/firedrake.git#[check,vtk,netgen,slepc]"
        
    # # Install slepc4py
    # export SLEPC_SRC_DIR="$(find $PETSC_DIR/$PETSC_ARCH/externalpackages \
    #                     -maxdepth 1 -name '*slepc*')"
    # pip install -vvvv --no-build-isolation --no-binary mpi4py,randomgen,islpy,numpy \
    #         --no-deps --ignore-installed $SLEPC_SRC_DIR/src/binding/slepc4py

    # # Install pkgs with option editable
    # git clone <fiat url>
    # pip install --editable ./fiat
}

curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-configure

# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# brew install $(python3 firedrake-configure --show-system-packages)
# brew install autoconf automake boost cmake gcc libtool make ninja openblas openmpi pkg-config python hwloc hdf5-mpi

# clone petsc
if [[ ! -d "./petsc" ]]; then
    git clone --depth 1 --branch $(python3 ./firedrake-configure --show-petsc-version) https://gitlab.com/petsc/petsc.git
fi
PETSC_DIR=$(readlink -f .)/petsc

echo "=== Install Real Version ==="
ARCH=default
build_petsc $ARCH $PETSC_EXTRA_REAL
install_firedrake $ARCH
echo "=== Install Real Version Done ==="

echo -e "\n\n"
echo "=== Install Complex Version ==="
ARCH=complex
build_petsc $ARCH $PETSC_EXTRA_COMPLEX
install_firedrake $ARCH
echo "=== Install Complex Version Done ==="
