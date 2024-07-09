1. Create directory to install firedrake

    yang@606d6cb4b0c4:~$ mkdir -p $HOME/opt/firedrake
    yang@606d6cb4b0c4:~$ cd $HOME/opt/firedrake

2. Download the install script

    yang@606d6cb4b0c4:~/opt/firedrake$ curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 94330  100 94330    0     0   142k      0 --:--:-- --:--:-- --:--:--  142k

3. Create a shell scipt `install.sh` with install command

    yang@606d6cb4b0c4:~/opt/firedrake$ cat install.sh
    PETSC_CONFIGURE_OPTIONS=" \
        --download-fftw --download-mmg --download-p4est --download-parmmg \
        --download-triangle --download-tetgen --download-ctetgen \
        --download-hpddm --download-libpng \
        --download-pragmatic --download-eigen" \
    python3 firedrake-install --disable-ssh \
        --documentation-dependencies \
        --netgen --slepc \
        --venv-name $HOME/opt/firedrake/firedrake-real-int32

4. Create the dir to store pkgs for petsc

    yang@606d6cb4b0c4:~/opt/firedrake$ mkdir $HOME/opt/firedrake/pkgs

5. Run the install script first time

   This process may take 1-2 minutes before raise exception as expected.
   Some output is omit as its too long.

    yang@606d6cb4b0c4:~/opt/firedrake$ bash install.sh
    Running firedrake-install --disable-ssh --documentation-dependencies --netgen --slepc --venv-name /home/yang/opt/firedrake/firedrake-real-int32
    Checking for presence of package build-essential...
      installed.
    Checking for presence of package autoconf...
      installed.
    Checking for presence of package automake...
      installed.
    ...
    ...
    Successfully cloned repository loopy.
    Checking out branch main
    Successfully checked out branch main
    Updating submodules.
    Successfully updated submodules.
    Installing petsc/
    Depending on your platform, PETSc may take an hour or more to build!
    Traceback (most recent call last):
      File "/home/yang/opt/firedrake/firedrake-install", line 1865, in <module>
        install("petsc/")
      File "/home/yang/opt/firedrake/firedrake-install", line 1055, in install
        build_and_install_petsc()
      File "/home/yang/opt/firedrake/firedrake-install", line 1167, in build_and_install_petsc
        check_call([python, "./configure", "PETSC_DIR={}".format(petsc_dir), "PETSC_ARCH={}".format(petsc_arch)] + petsc_options)
      File "/home/yang/opt/firedrake/firedrake-install", line 668, in check_call
        log.debug(subprocess.check_output(arguments, stderr=subprocess.STDOUT, env=os.environ).decode())
      File "/usr/lib/python3.10/subprocess.py", line 421, in check_output
        return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
      File "/usr/lib/python3.10/subprocess.py", line 526, in run
        raise CalledProcessError(retcode, process.args,
    subprocess.CalledProcessError: Command '['/home/yang/opt/firedrake/firedrake-real-int32/bin/python', './configure', 'PETSC_DIR=/home/yang/opt/firedrake/firedrake-real-int32/src/petsc', 'PETSC_ARCH=default', '--download-pnetcdf', '--download-metis', '--download-tetgen', '--download-chaco', '--with-c2html=0', '--download-p4est', '--download-pragmatic', '--download-mmg', '--with-debugging=0', '--download-hdf5', '--with-fortran-bindings=0', '--download-scalapack', '--with-packages-download-dir=/home/yang/opt/firedrake/pkgs', '--download-hpddm', '--download-eigen', '--download-bison', '--with-shared-libraries=1', '--download-mumps', '--with-zlib', '--download-hwloc', '--download-libpng', '--download-triangle', '--download-ptscotch', '--download-mpich', '--download-superlu_dist', '--download-cmake', '--download-fftw', '--download-hypre', '--download-suitesparse', '--download-ctetgen', '--download-pastix', '--download-netcdf', '--download-parmmg']' returned non-zero exit status 10.


    Install log saved in /home/yang/opt/firedrake/firedrake-install.log


6. Obtain the packages to download from the log file `firedrake-install.log`

    yang@606d6cb4b0c4:~/opt/firedrake$ sed -n '/Download the following/,/Then run/p' firedrake-install.log > pkgs_info.txt
    yang@606d6cb4b0c4:~/opt/firedrake$ cat pkgs_info.txt
    2024-02-01 10:37:28,382 DEBUG  Download the following packages to /home/yang/opt/firedrake/pkgs

    mpich ['https://github.com/pmodels/mpich/releases/download/v4.1.2/mpich-4.1.2.tar.gz', 'https://www.mpich.org/static/downloads/4.1.2/mpich-4.1.2.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/mpich-4.1.2.tar.gz']
    hwloc ['https://download.open-mpi.org/release/hwloc/v2.10/hwloc-2.10.0.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/hwloc-2.10.0.tar.gz']
    cmake ['https://github.com/Kitware/CMake/releases/download/v3.28.1/cmake-3.28.1.tar.gz', 'https://gitlab.kitware.com/cmake/cmake/-/archive/v3.28.1/cmake-v3.28.1.tar.gz']
    hdf5 ['https://web.cels.anl.gov/projects/petsc/download/externalpackages/hdf5-1.14.3-p1.tar.bz2']
    netcdf ['https://web.cels.anl.gov/projects/petsc/download/externalpackages/netcdf-c-4.9.2-p1.tar.gz']
    pnetcdf ['git clone https://github.com/parallel-netcdf/pnetcdf', 'https://github.com/parallel-netcdf/pnetcdf/archive/checkpoint.1.12.3.tar.gz', 'https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/pnetcdf-1.12.3.tar.gz']
    hypre ['git clone https://github.com/hypre-space/hypre', 'https://github.com/hypre-space/hypre/archive/v2.30.0.tar.gz']
    chaco ['git clone https://bitbucket.org/petsc/pkg-chaco.git', 'https://bitbucket.org/petsc/pkg-chaco/get/v2.2-p4.tar.gz']
    metis ['git clone https://bitbucket.org/petsc/pkg-metis.git', 'https://bitbucket.org/petsc/pkg-metis/get/v5.1.0-p11.tar.gz']
    suitesparse ['git clone https://github.com/DrTimothyAldenDavis/SuiteSparse', 'https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v7.6.0.tar.gz']
    eigen ['git clone https://gitlab.com/libeigen/eigen.git', 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/eigen-3.4.0.tar.gz']
    ptscotch ['git clone https://gitlab.inria.fr/scotch/scotch.git', 'https://gitlab.inria.fr/scotch/scotch/-/archive/v7.0.3/scotch-v7.0.3.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/scotch-v7.0.3.tar.gz']
    bison ['https://ftp.gnu.org/gnu/bison/bison-3.8.2.tar.gz', 'http://mirrors.kernel.org/gnu/bison/bison-3.8.2.tar.gz']
    mumps ['https://graal.ens-lyon.fr/MUMPS/MUMPS_5.6.2.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/MUMPS_5.6.2.tar.gz']
    scalapack ['git clone https://github.com/Reference-ScaLAPACK/scalapack', 'https://github.com/Reference-ScaLAPACK/scalapack/archive/5bad7487f496c811192334640ce4d3fc5f88144b.tar.gz']
    pastix ['https://web.cels.anl.gov/projects/petsc/download/externalpackages/pastix_5.2.3.tar.bz2']
    superlu_dist ['git clone https://github.com/xiaoyeli/superlu_dist', 'https://github.com/xiaoyeli/superlu_dist/archive/v8.2.1.tar.gz']
    triangle ['git clone https://bitbucket.org/petsc/pkg-triangle', 'https://bitbucket.org/petsc/pkg-triangle/get/v1.3-p2.tar.gz']
    ctetgen ['git clone https://bitbucket.org/petsc/ctetgen', 'https://bitbucket.org/petsc/ctetgen/get/ctetgen-0.11.tar.gz']
    fftw ['https://www.fftw.org/fftw-3.3.10.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/fftw-3.3.10.tar.gz']
    hpddm ['git clone https://github.com/hpddm/hpddm', 'https://github.com/hpddm/hpddm/archive/201eecd26177f88d7bb6287251877d8013fb64d2.tar.gz']
    libpng ['https://sourceforge.net/projects/libpng/files/libpng16/1.6.40/libpng-1.6.40.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/libpng-1.6.40.tar.gz']
    mmg ['git clone https://github.com/MmgTools/mmg.git', 'https://github.com/MmgTools/mmg/archive/cc54b4174871212cd32595c5dca732d40e01e90b.tar.gz']
    p4est ['git clone https://github.com/cburstedde/p4est', 'https://github.com/cburstedde/p4est/archive/56b58bd7a5462ef85e136cea9fd9ee6bf9558e71.tar.gz']
    parmmg ['git clone https://github.com/MmgTools/ParMmg.git', 'https://github.com/MmgTools/ParMmg/archive/eaebcbe53488f11ce68325502c9f760ed250e0b3.tar.gz']
    pragmatic ['git clone https://github.com/meshadaptation/pragmatic.git', 'https://github.com/meshadaptation/pragmatic/archive/ef941eddf50a6de307a5d6b54b5d44504dd3ce89.tar.gz']
    tetgen ['http://www.tetgen.org/1.5/src/tetgen1.6.0.tar.gz', 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/tetgen1.6.0.tar.gz']

    Then run the script again

7. Download the pkgs
7.1 Download the script to download pkgs
    yang@606d6cb4b0c4:~/opt/firedrake$ curl -L -O https://zzyang.net/firedrake-notes/_downloads/15cd80c553ab2d83b1094de3a5a5aa5c/download_petsc_pkgs.py
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100  2930  100  2930    0     0   4494      0 --:--:-- --:--:-- --:--:--  4493

7.2 Run the download script

    This process may take several minutes depending on the network

    yang@606d6cb4b0c4:~/opt/firedrake$ python3 download_petsc_pkgs.py -d $HOME/opt/firedrake/pkgs pkgs_info.txt
    Downloading 27 packages to /home/yang/opt/firedrake/pkgs...
    [01/27] Download mpich...
      [1/3] curl -L -O https://github.com/pmodels/mpich/releases/download/v4.1.2/mpich-4.1.2.tar.gz
    [02/27] Download hwloc...
      [1/2] curl -L -O https://download.open-mpi.org/release/hwloc/v2.10/hwloc-2.10.0.tar.gz
    ......
    [26/27] Download pragmatic...
      [1/2] git clone https://github.com/meshadaptation/pragmatic.git
    [27/27] Download tetgen...
      [1/2] curl -L -O http://www.tetgen.org/1.5/src/tetgen1.6.0.tar.gz

    All packages downloaded successfully.

8. Clean the installation of the first time

    yang@606d6cb4b0c4:~/opt/firedrake$ rm -rf firedrake-real-int32

9. Run the install script for the second time

    yang@606d6cb4b0c4:~/opt/firedrake$ bash install.sh

10. Activate the env to use firedrake
