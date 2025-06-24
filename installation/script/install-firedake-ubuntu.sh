#!/bin/bash

set -e

ARCH=${1:-real-cuda} # can be real, real-int64, real-cuda, complex, complex-int64, complex-cuda
if [ $# -ge 1 ]; then
  shift
fi

ssh -v -n -N -D 5000 -o ExitOnForwardFailure=yes proxy-server >ssh-proxy.log 2>&1 &
SSH_PID=$!

function set_proxy() {
    MY_PROXY="socks5h://127.0.0.1:5000"
    export http_proxy=$MY_PROXY
    export https_proxy=$MY_PROXY
}

function unset_proxy() {
    unset http_proxy https_proxy
}


for i in {1..20}; do
    if ! ps -p $SSH_PID > /dev/null; then
        echo "SSH failed to start. See $LOGFILE for details."
        exit 1
    fi
    sleep 0.1
done

echo "Started SSH tunnel with PID $SSH_PID"

trap "echo 'Killing SSH tunnel...'; kill $SSH_PID" EXIT

PETSC_PKG_DOWNLOAD_DIR="--with-packages-download-dir=/opt/firedrake/petsc-pkgs"

PETSC_OPTIONS_BASE="
  --download-slepc \
"

PETSC_INT64="--with-64-bit-indices"
PETSC_INT32="--download-parmmg"

PETSC_EXTRA_COMMON=" \
    --download-hpddm --download-libpng \
    --download-ctetgen --download-tetgen --download-triangle \
    --download-mmg --download-p4est"

PETSC_EXTRA_REAL="--download-eigen --download-hypre"
PETSC_EXTRA_REAL_INT32="--download-pragmatic"
PETSC_EXTRA_COMPLEX="--with-scalar-type=complex"
PETSC_EXTRA_CUDA="--with-cuda-dir=/usr/local/cuda --download-kokkos --download-kokkos-kernels"

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
    set_proxy
    time ./configure PETSC_ARCH=$PETSC_ARCH \
                --with-c2html=0 \
                --with-debugging=0 \
                --with-fortran-bindings=0 \
                --with-shared-libraries=1 \
                --with-strict-petscerrorcode \
                --COPTFLAGS="-O3 -march=native -mtune=native \
                             -I/usr/include/hdf5/openmpi \
                             -I/usr/include/scotch \
                             -I/usr/include/superlu \
                             -I/usr/include/superlu-dist \
                             -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi" \
                --CXXOPTFLAGS="-O3 -march=native -mtune=native \
                             -I/usr/include/hdf5/openmpi \
                             -I/usr/include/scotch \
                             -I/usr/include/superlu \
                             -I/usr/include/superlu-dist \
                             -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi" \
                --FOPTFLAGS="-O3 -march=native -mtune=native \
                             -I/usr/include/hdf5/openmpi \
                             -I/usr/include/scotch \
                             -I/usr/include/superlu \
                             -I/usr/include/superlu-dist \
                             -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi" \
                --with-bison \
                --with-fftw \
                --with-hdf5 \
                --with-hwloc \
                --with-metis \
                --with-mumps \
                --with-netcdf \
                --with-pnetcdf \
                --with-ptscotch \
                --with-scalapack-lib=-lscalapack-openmpi \
                --with-suitesparse \
                --with-superlu_dist \
                --with-zlib \
                "$@"
    time make PETSC_ARCH=$PETSC_ARCH all
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
    
    if [[ "$HPC_LOGIN_NODE" != "" ]]; then
        ssh $HPC_LOGIN_NODE "source $(pwd)/$ENV_NAME/bin/activate && pip install pysocks"
    else
        unset_proxy
        pip install pysocks
    fi

    echo ENV: CC=mpicc CXX=mpicxx
    echo ENV: PETSC_DIR=$PETSC_DIR
    echo ENV: PETSC_ARCH=$PETSC_ARCH
    echo ENV: SLEPC_DIR=$SLEPC_DIR 
    echo ENV: HDF5_MPI=ON # HDF5_DIR=/opt/homebrew

    set_proxy
    CC=mpicc CXX=mpicxx \
        PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH SLEPC_DIR=$SLEPC_DIR \
        HDF5_MPI=ON `#HDF5_DIR=$PETSC_DIR/$PETSC_ARCH` \
        time pip install -v --no-binary h5py --src $ENV_NAME/src \
            --editable "git+https://github.com/firedrakeproject/firedrake.git@release#egg=firedrake[check,vtk,netgen,slepc]"
    PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH SLEPC_DIR=$SLEPC_DIR \
    time firedrake-check
}

set_proxy
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-configure

# if [[ ! -d "./petsc" ]]; then
#     git clone --depth 1 --branch $(python3 ./firedrake-configure --show-petsc-version) https://gitlab.com/petsc/petsc.git
# fi

# clone petsc
PETSC_DIR=$(readlink -f .)/petsc

PETSC_OPTIONS="$PETSC_PKG_DOWNLOAD_DIR $PETSC_OPTIONS_BASE $PETSC_EXTRA_COMMON"

if [[ "$ARCH" == *default* || "$ARCH" == *real* ]]; then
    PETSC_OPTIONS="$PETSC_OPTIONS $PETSC_EXTRA_REAL"
elif [[ "$ARCH" == *complex* ]]; then
    PETSC_OPTIONS="$PETSC_OPTIONS $PETSC_EXTRA_COMPLEX"
fi

if [[ "$ARCH" == *int64* ]]; then
    PETSC_OPTIONS="$PETSC_OPTIONS $PETSC_INT64"
else # int32
    if [[ "$ARCH" == *default* || "$ARCH" == *real* ]]; then
        # real int32
        PETSC_OPTIONS="$PETSC_OPTIONS $PETSC_EXTRA_REAL_INT32"
    fi
fi

if [[ "$ARCH" == *cuda* ]]; then
    PETSC_OPTIONS="$PETSC_OPTIONS $PETSC_EXTRA_CUDA"
fi


echo "=== Install firedrake $ARCH Version ==="
echo `date` run build_petsc $ARCH $PETSC_OPTIONS "$@"
build_petsc $ARCH $PETSC_OPTIONS "$@"
echo `date` run install_firedrake $ARCH
install_firedrake $ARCH
echo "=== Install firedrake $ARCH Version Done ==="