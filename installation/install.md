# [Installation of Firedake](https://www.firedrakeproject.org/download.html)


<!--
> {sub-ref}`today` | {sub-ref}`wordcount-words` words | {sub-ref}`wordcount-minutes` min read

{attribution="Hamlet act 4, Scene 5"}
> We know what we are, but know not what we may be.
-->

To install firedrake, the computer should have access to the Internet.
Otherwise, please refer to {doc}`install_without_internet_access`.

## Ubuntu

### Install with default options

The easiest way to intall firedrake is to download the installation script `firedrake-install` and run it using Python.
This method will intall the real number version by defaults.

1. Download the installation script

   ```bash
   curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   ```

2. Install by running the script

   ```bash
   python3 firedrake-install
   ```

   If you need to know more installation options, please refer to the help documentation.

   ```bash
   python3 firedrake-install -h
   ```

````{note}
Sometimes there may be issues with accessing the `pip` source during installation, and error messages similar to the following may appear:

  ```console
  Starting new HTTPS connection (6): pypi.org:443
  Could not fetch URL https://pypi.org/simple/pulp/: connection error: HTTPSConnectionPool(host='pypi.org', port=443): Max retries exceeded with url: /simple/pulp/ (Caused by NewConnectionError('<pip._vendor.urllib3.connection.HTTPSConnection object at 0x7f43dce52bc0>: Failed to establish a new connection: [Errno 101] Network is unreachable')) - skipping
  Skipping link: not a file: https://pypi.org/simple/pulp/
  Given no hashes to check 0 links for project 'pulp': discarding no candidates
  ERROR: Could not find a version that satisfies the requirement PuLP (from versions: none)
  ERROR: No matching distribution found for PuLP
  ```

This can be fixed by setting the source of `pip`, such as changing it to the source of USTC:

```bash
mkdir -p $HOME/.pip && \
cat > $HOME/.pip/pip.conf <<EOF
[global]
index-url = https://pypi.mirrors.ustc.edu.cn/simple
[install]
trusted-host=pypi.mirrors.ustc.edu.cn
EOF
```
````

````{note}
After the installation, we can install some useful python packages:
```bash
pip install gmsh, meshio
```
````

### Test and run examples

#### Test firedrake

```bash
source <path-to-firedrake-env>/bin/activate
cd $VIRTUAL_ENV/src/firedrake
pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
```

#### Run examples

Create a python file `example.py` with the following content:

```python
from firedrake import *

N = 8
mesh = RectangleMesh(nx=N, ny=N, Lx=1, Ly=1)
x, y = SpatialCoordinate(mesh)
f = sin(pi*x)*sin(pi*y)
g = Constant(0)

V = FunctionSpace(mesh, 'CG', degree=1)

u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx

bc = DirichletBC(V, g=g, sub_domain='on_boundary')

u_h = Function(V, name='u_h')
solve(a == L, u_h, bcs=bc)
```

Run in serial:
```bash
python example.py
```

Run in parallel:
```bash
mpiexec -n 4 python example.py
```

### Install with custom options

#### Install `real-int32` and/or `real-int32-debug`

1. Download the installation script

   ```bash
   curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   ```

2. Enable PETSc's debug option (optional)

   This option is turned off in firedrake by default. If you need to turn it on, you can use the following command to enable it.

   ```bash
   DEBUG='-debug'
   sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' firedrake-install
   ```

3. Install

   ```bash
   PETSC_CONFIGURE_OPTIONS=" \
       --download-fftw --download-mmg --download-p4est --download-parmmg \
       --download-triangle --download-tetgen --download-ctetgen \
       --download-hpddm --download-libpng \
       --download-pragmatic --download-eigen" \
   python3 firedrake-install --disable-ssh \
       --documentation-dependencies \
       --netgen --slepc \
       --venv-name $HOME/opt/firedrake/firedrake-real-int32$DEBUG
   ```

   ```{note}
   If you would like to install int64 version, please use `--petsc-int-type int64` instead.

   `pragmatic` cannot be used with `int64`
   ```


#### Install `complex-int32`, and/or `complex-int32-debug`

1. Download the installation script

   ```bash
   curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   ```

2. Enable PETSc's debug option (optional)

   ```bash
   DEBUG='-debug'
   sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' firedrake-install
   ```

3. Install

   ```bash
   PETSC_CONFIGURE_OPTIONS=" \
       --download-fftw --download-mmg --download-p4est --download-parmmg \
       --download-triangle --download-tetgen --download-ctetgen \
       --download-hpddm --download-libpng" \
   python3 firedrake-install --disable-ssh \
       --documentation-dependencies  \
       --netgen --slepc \
       --petsc-int-type int32 --complex \
       --venv-name $HOME/opt/firedrake/firedrake-complex-int32$DEBUG
    ```

#### Install with MKL

If you want to use solvers provided by MKL, you can follow the following steps to install firedrake with MKL

```{warning}
The following steps are not tested for the recent version of firedrake.
```

1. Install mkl

   1. Add repo of mkl

      ```bash
      wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
      | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

      echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] \
      https://apt.repos.intel.com/oneapi all main" \
      | sudo tee /etc/apt/sources.list.d/oneAPI.list

      sudo apt update
      ```

   2. Install libs and headers of MKL

      ```bash
      # sudo apt install intel-basekit
      sudo apt install intel-oneapi-mkl
      sudo apt install intel-oneapi-mkl-devel
      ```

3. Download the installation script and enable the debug option if necessary

   1. Download the script as before

      ```bash
      curl -O \
      https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
      ```

   2. Enable the debug option and patch the script for `mkl` and `hpddm`

      ```bash
      sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' -e 's/\({0}\/lib\)/\1\/intel64/g' \
          -e 's/\(.*\)\(--C\)\(FLAGS=-I{}\/include\)\(.*\)/\1\2\3\4\n\1\2XX\3\4/' \
          firedrake-install
      ```

      - `'s/\(--with-debugging=\)0/\11/g'` for petsc debug
      - `'s/\({0}\/lib\)/\1\/intel64/g'` for `mkl` lib
      - `'s/\(.*\)\(--C\)\(FLAGS=-I{}\/include\)\(.*\)/\1\2\3\4\n\1\2XX\3\4/'` for `hpddm` with `mkl`

4. Install Firedrake `real-int32`

   ```bash
   time PETSC_CONFIGURE_OPTIONS="--download-fftw --download-mmg \
           --download-p4est --download-parmmg --download-triangle \
           --download-tetgen --download-ctetgen --download-hpddm --download-libpng \
           --download-slepc  --download-pragmatic --download-eigen \
           --with-mkl_pardiso-dir=/opt/intel/oneapi/mkl/latest \
           --with-mkl_cpardiso-dir=/opt/intel/oneapi/mkl/latest" \
   python3 firedrake-install --disable-ssh --documentation-dependencies \
       --with-blas=/opt/intel/oneapi/mkl/latest \
       --venv-name $HOME/opt/firedrake/firedrake-real-int32-mkl-debug
   ```

5. Fix the error on `mkl_cpardiso`

   If you run the following test, there will be an error:

   ```console
   $ cd $(dirname `which python`)/../src/petsc/src/binding/petsc4py/demo/kspsolve
   $ python test_mat_ksp.py -pc_type lu -pc_factor_mat_solver_type mkl_cpardiso -ksp_view
   Intel MKL FATAL ERROR: Cannot load symbol MKLMPI_Get_wrappers.
   ```

   1. Patch `petsc4py`

      Patch `petsc4py` by the following file

      ```console
      $ git diff
      diff --git a/src/binding/petsc4py/conf/confpetsc.py b/src/binding/petsc4py/conf/confpetsc.py
      index 5801b146ff..b00fab2d32 100644
      --- a/src/binding/petsc4py/conf/confpetsc.py
      +++ b/src/binding/petsc4py/conf/confpetsc.py
      @@ -319,6 +319,11 @@ class PetscConfig:
               self._configure_ext(extension, petsc_inc, preppend=True)
               self._configure_ext(extension, petsc_lib)

      +        blas_lib = flaglist(self['BLASLAPACK_LIB'])
      +        blas_inc = flaglist(self['BLASLAPACK_INCLUDE'])
      +        self._configure_ext(extension, blas_inc, preppend=True)
      +        self._configure_ext(extension, blas_lib)
      +
           def configure_compiler(self, compiler):
               if compiler.compiler_type != 'unix': return
               getenv = os.environ.get
      ```

      ````{note}
      The value of BLASLAPACK_LIB is

      ```bash
      LASLAPACK_LIB="-Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib/intel64 \
          -L/opt/intel/oneapi/mkl/latest/lib/intel64 \
          -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread \
          -lmkl_blacs_intelmpi_lp64 -lgomp -ldl -lpthread"
      ```

      The other way to fix this is modifing the file `firedrake-install` by adding the following content to `blas["LDFLAGS"]`

      ```python
      blas["LDFLAGS"] = "-Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib/intel64 \
          -L/opt/intel/oneapi/mkl/latest/lib/intel64 \
          -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread \
          -lmkl_blacs_intelmpi_lp64 -lgomp -ldl -lpthread"
      ```
      ````

   2. Recompile and install petsc4py (in the activated Firedrake environment)

      ```bash
      export PETSC_DIR=$(readlink -f $(dirname `which python`)/../src/petsc)
      export PETSC_ARCH=default
      cd $PETSC_DIR/src/binding/petsc4py
      make clean
      python -m pip install --no-build-isolation --no-binary mpi4py,randomgen,islpy,numpy \
          --no-deps -vvv --ignore-installed .
      ```

6. Install `slepc4py`

   ```bash
   export SLEPC_DIR="$(find $PETSC_DIR/$PETSC_ARCH/externalpackages \
                       -maxdepth 1 -name '*slepc*')"
   python -m pip install --no-build-isolation --no-binary mpi4py,randomgen,islpy,numpy \
       --no-deps -vvv --ignore-installed $SLEPC_DIR/src/binding/slepc4py
   ```

The complex version `complex-int32` with debug can be installed using the following command.

```bash
PETSC_CONFIGURE_OPTIONS=" \
    --download-fftw --download-mmg --download-pragmatic --download-eigen \
    --download-p4est --download-parmmg --download-triangle \
    --download-tetgen --download-ctetgen --download-hpddm --download-libpng \
    --download-slepc --download-scalapack --download-mumps \
    --with-mkl_pardiso-dir=/opt/intel/oneapi/mkl/latest \
    --with-mkl_cpardiso-dir=/opt/intel/oneapi/mkl/latest" \
python3 firedrake-install --disable-ssh \
    --documentation-dependencies  \
    --with-blas=/opt/intel/oneapi/mkl/latest --complex \
    --venv-name firedrake/complex-int32-mkl-debug
```

If you encounter the same error of solver `mkl_cpardiso`, you can fix it by using the same method as before.

#### Install with petsc packages pre-downloaded into a directory

Sometimes, some of the packages that petsc depends on cannot be downloaded automatically.
We can add the option

```bash
--with-packages-download-dir=<path/to/petsc/packages>
```

to command `configure` of petsc to obtain the list of required packages.
Then download these packages manually and put them into the path.
Afterwards, configure it again with the above option.

The python script [`download_petsc_pkgs.py`](./script/download_petsc_pkgs.py) can be used to download the packages.
Save the output of petsc configure to a file, for example [`pkgs_info.txt`](./script/pkgs_info.txt).
Then run the following command to download packages:

```bash
python3 download_petsc_pkgs.py -d <path/to/petsc/packages> pkgs_info.txt
```

Here is an example illustrating how to install firedrake with the above option.

1. Create directory to install firedrake

    ```bash
    yang@606d6cb4b0c4:~$ mkdir -p $HOME/opt/firedrake
    yang@606d6cb4b0c4:~$ cd $HOME/opt/firedrake
    ```

2. Download the install script

    ```bash
    yang@606d6cb4b0c4:~/opt/firedrake$ curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 94330  100 94330    0     0   142k      0 --:--:-- --:--:-- --:--:--  142k
    ```

3. Create a shell scipt `install.sh` with install command

    ```bash
    yang@606d6cb4b0c4:~/opt/firedrake$ cat install.sh
    PETSC_CONFIGURE_OPTIONS=" \
        --with-packages-download-dir=$HOME/opt/firedrake/pkgs
        --download-fftw --download-mmg --download-p4est --download-parmmg \
        --download-triangle --download-tetgen --download-ctetgen \
        --download-hpddm --download-libpng \
        --download-pragmatic --download-eigen" \
    python3 firedrake-install --disable-ssh \
        --documentation-dependencies \
        --netgen --slepc \
        --venv-name $HOME/opt/firedrake/firedrake-real-int32
    ```

4. Create the dir to store pkgs for petsc

    ```bash
    yang@606d6cb4b0c4:~/opt/firedrake$ mkdir $HOME/opt/firedrake/pkgs
    ```

5. Run the install script first time

   This process may take 1-2 minutes before raising exception as expected.
   Some output is omit as its too long.

    ```bash
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
    ```


6. Obtain the packages to download from the log file `firedrake-install.log`

    ```bash
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
    ```

7. Download the pkgs

   1. Download the script to download pkgs

      ```bash
      yang@606d6cb4b0c4:~/opt/firedrake$ curl -L -O https://zzyang.net/firedrake-notes/_downloads/15cd80c553ab2d83b1094de3a5a5aa5c/download_petsc_pkgs.py
         % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                       Dload  Upload   Total   Spent    Left  Speed
      100  2930  100  2930    0     0   4494      0 --:--:-- --:--:-- --:--:--  4493
      ```

   2. Run the download script

      This process may take several minutes depending on the network

      ```bash
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
      ```

8. Clean the installation of the first time

    ```bash
    yang@606d6cb4b0c4:~/opt/firedrake$ rm -rf firedrake-real-int32
    ```

9. Run the install script for the second time

   This may take half to one hour to finish depending on the performance of the computer and the network.

    ```bash
    yang@606d6cb4b0c4:~/opt/firedrake$ bash install.sh
    ```

10. Activate the env to use firedrake according to the output of the last step

    ```bash
    yang@606d6cb4b0c4:~/opt/firedrake$ source $HOME/opt/firedrake/firedrake-real-int32/bin/activate
    ```

### Some notes on PETSc

1. Enable PETSc with X

   1. Install `libx11-dev`
   
      ```bash
      sudo apt install libx11-dev
      ```
   
   2. Enable `--with-x=1` option

      This option is turned off in firedrake by default.
      If you need to turn it on, you can use the following command to enable it before run the install script.

      ```bash
      sed -i.bak -e 's/\(--with-x=\)0/\11/g' firedrake-install
      ```

2. Add `bin` path of petsc to `PATH`
   
   The `bin` directory of petsc provides some useful tools, but it is not added to the `PATH` by default.
   We can define command `add-petsc-bin`.

   ```bash
   alias add-petsc-bin='export \
       PATH=$PATH:$(dirname $(which python))/../src/petsc/lib/petsc/bin:$(\
       dirname $(which python))/../src/petsc/default/bin'
   ```

   Executing it in activated firedrake env will add the petsc/bin to `PATH`.


### Using firedrake in Jupyter-lab

1. Install `jupyterlab`

    ```bash
    python3 -m pip install jupyterlab
    ```

    Maybe you need add `$HOME/.local/bin` to environment variable `PATH`:

    ```bash
    export PATH=$PATH:$HOME/.local/bin
    ```

2. Configure jupyterlab

    Generate config file for customizing jupyter lab:

    ```bash
    jupyter lab --generate-config
    ```

    Set password for jupyter lab:

    ```bash
    jupyter lab password  #
    ```

3. Configure Browser

    In wsl-ubuntu, configure the browser like this:

    ```bash
    export BROWSER="/path/to/chrome/or/firefox"
    ```

    An example of chrome:

    ```bash
    export BROWSER='/mnt/c/Program Files/Google/Chrome/Application/chrome.exe'
    ```

    Now, you can type `jupyter-lab` to start jupyter. You will see jupyter in browser.

4. Configure kernels

    1. Activate env:

        ```console
        $ source /your/env/path/activate
        ```

    2. Add kernels:

        ```console
        (your-venv)$ ipython kernel install --name "local-venv" --user
        ```

        The output should looks like:
        ```console
        Installed kernelspec local-venv in /path/to/kernels/local-venv
        ```
        where the path `/path/to/kernels/local-venv` should have different name in your enviornment.

        Also, we can find the kernelspec path by run `jupyter kernelspec list`, which will output all avaialbe kernels with its path.

        Now we need check the python path in `kernel.json`, which is located in the above path.
        Make sure it is the python in your env. Otherwise, correct it.

    3. Add environment variables to `kernel.json`:[^kernels]

       [^kernels]: For more about jupyter kernels: https://jupyter-client.readthedocs.io/en/stable/kernels.html

       An exmaple of `kernel.json`:

        ```json
        {
         "argv": [
          "/home/yzz/firedrake/real-int32-debug/bin/python",
          "-m",
          "ipykernel_launcher",
          "-f",
          "{connection_file}"
         ],
         "env": {
          "OMP_NUM_THREADS": "1",
          "PATH": "/home/yzz/firedrake/real-int32-debug/bin:${PATH}"
         },
         "display_name": "firedrake-real-int32",
         "language": "python",
         "metadata": {
          "debugger": true
         }
        }
        ```

### Update firedrake

Generally, you can simply run `firedrake-update` in the activated environment to update firedrake.

If you want to rebuild PETSc (i.e., using the `--rebuild` option)
and you have used `PETSC_CONFIGURE_OPTIONS` and `--with-blas` during the installation,
you also need to use these two options when updating.

In addition, if you installed MKL using the aforementioned method,
you may need to modify the firedrake-update script.

The example for `firedrake/complex-int32-mkl-debug` is as follows:

1. Modify `firedrake-update`

    ```bash
    sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' -e 's/\({0}\/lib\)/\1\/intel64/g' \
        -e 's/\(.*\)\(--C\)\(FLAGS=-I{}\/include\)\(.*\)/\1\2\3\4\n\1\2XX\3\4/' \
        firedrake-update
    ```

2. Update

    ```bash
    PETSC_CONFIGURE_OPTIONS=" \
        --download-fftw --download-mmg --download-pragmatic --download-eigen \
        --download-p4est --download-parmmg --download-triangle \
        --download-tetgen --download-ctetgen --download-hpddm --download-libpng \
        --download-slepc --download-scalapack --download-mumps \
        --with-mkl_pardiso-dir=/opt/intel/oneapi/mkl/latest \
        --with-mkl_cpardiso-dir=/opt/intel/oneapi/mkl/latest" \
    firedrake-update --rebuild --no-update-script --with-blas=/opt/intel/oneapi/mkl/latest
    ```

## Windows

Install WSL (Windows Subsystem for Linux) on Windows (the system installed is `Ubuntu` by default) and then install Firedrake as before.

### Install WSL

To install WSL, please see https://docs.microsoft.com/zh-cn/windows/wsl/install.

### Install Firedrake

Follow the installation method for Ubuntu.

<!--
### Mount network locations in wsl-ubuntu

1. Map the network locations to a driver, for example to `Y:`

2. Edit `/etc/fstab` as follows (`sudo vi /etc/fstab`):

    ```bash
    $ cat /etc/fstab
    # UNCONFIGURED FSTAB FOR BASE SYSTEM
    Y: /mnt/y drvfs auto,users,dev,exec,rw,async,relatime,uid=1000,gid=1000 0 0
    ```

3. Mount:

    ```bash
    sudo mount -a
    ```
-->

## MacOS

First, install Homebrew[^brew], and then use Homebrew to install python3. After that, install Firedrake directly, similar to Ubuntu.

Please install `pkgconf` before start the installation if you add package `p4est` for PETSc.

[^brew]: Homebrew https://brew.sh/

1. Install Homebrew

   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

2. Install python3

   When installing python3, you can choose the version you want to install. For example, you can install python3.11 as follows.

   ```bash
   brew install python3
   brew install m4
   ```

   ````{note}
   If you using Xcode command line tool other then Xcode. You should install m4 using brew or patch the system with m4.
   ```
   xcode-select -p                      # this will print the path to command line tools
   cd <the-path-to-command-line-tools>  # should be /Library/Developer/CommandLineTools
   cd usr/bin
   ln -s gm4 m4                         # need password for the system
   ```
   ````

3. Install firedrake

   Now, you can follow the installation method for Ubuntu.

   ```bash
   curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   $(brew --prefix)/bin/python3 firedrake-install
   ```


## Linux Server

If the server cannot access the network, please refer {doc}`install_without_internet_access`.
The Firedrake team provides a way[^firedrake-spack] to install Firedrake based on [Spack](https://github.com/spack/spack)[^spack], a package manager for HPC.

[^spack]: spack https://github.com/spack/spack.
[^firedrake-spack]: firedrake-spack https://github.com/firedrakeproject/firedrake-spack.


1. Download spack

   ```bash
   mkdir -p $HOME/opt
   cd $HOME/opt && \
   git clone -c feature.manyFiles=true https://github.com/lrtfm/spack.git && \
   pushd spack
   git checkout lrtfm/develop
   popd
   source $HOME/opt/spack/share/spack/setup-env.sh
   ```

   ````{note}
   Add the following command to the file $HOME/.bashrc to add shell support for spack.

   ```bash
   source $HOME/opt/spack/share/spack/setup-env.sh
   ```
   ````

   ````{note}
   On some workstations, the content of the /tmp directory may not have execution permissions. You need to change the spack build directory as follows.

   ```bash
   mkdir -p $HOME/.spack && \
   cat > $HOME/.spack/config.yaml <<EOF
   config:
     build_stage:
       - \$user_cache_path/stage
   EOF
   ```
   ````

2. Download firedrake-spack

   ```bash
   cd $HOME/opt && \
   git clone https://github.com/lrtfm/firedrake-spack.git && \
   pushd firedrake-spack && \
   git checkout lrtfm/air-gapped-install && \
   popd
   ```

   ```{note}
   The current version of petsc in firedrakeproject will break when using some compilers: https://lists.mcs.anl.gov/pipermail/petsc-users/2023-April/048482.html. Patch has been added branch `lrtfm/air-gapped-install` of `firedrake-spack`.
   ```

3. Create spack env and add packages

   + complex-int32

     1. Create spack env

        ```bash
        cd $HOME/opt && \
        FIREDRAKE_ENV_NAME=firedrake-complex-int32 && \
        spack env create -d $FIREDRAKE_ENV_NAME && \
        spack env activate -p $FIREDRAKE_ENV_NAME && \
        spack -e $SPACK_ENV config add concretizer:unify:true
        ```

     2. Add firedrake repo

        We add the firedrake repo to the created space env

        ```bash
        cd $HOME/opt && \
        spack repo add firedrake-spack
        ```

     3. Add packages

        ```bash
        spack add py-firedrake@develop%gcc +complex ^mpich ^openblas ^slepc+hpddm \
            ^petsc+libpng+libyaml+parmmg+mmg+hpddm+tetgen+valgrind \
            ^hypre+superlu-dist ^vtk@9.0.3 && \
        spack add gmsh py-meshio py-tqdm py-pyyaml py-memory-profiler
        ```

   + real-int32

     1. Create env

        ```bash
        cd $HOME/opt && \
        FIREDRAKE_ENV_NAME=firedrake-real-int32 && \
        spack env create -d $FIREDRAKE_ENV_NAME && \
        spack env activate -p $FIREDRAKE_ENV_NAME && \
        spack -e $SPACK_ENV config add concretizer:unify:true
        ```

     2. Add firedrake repo

        ```bash
        cd $HOME/opt && \
        spack repo add firedrake-spack
        ```

     3. Add packages

        Copy the follwing command into bash will raise error. Please copy line by line.

        ```bash
        spack add py-firedrake@develop%gcc ^mpich ^openblas ^slepc+hpddm \
            ^petsc+libpng+libyaml+parmmg+mmg+hpddm+tetgen+valgrind \
            ^hypre+superlu-dist ^vtk@9.0.3 && \
        spack add gmsh py-meshio py-tqdm py-pyyaml py-memory-profiler
        ```

   ```{note}
   Installing `vtk@8.x.x` and `vtk@9.2.2`(on some hosts) in spack will fail. We use `vtk@9.0.3`.(2023-04-30)
   ```

<!--
This is not needed. As mpich support slurm by default with hydary.

    __Remark 4__: "On a cluster that uses slurm to submit jobs, you need to add slurm as an external package, and then change `^mpich` to `^mpich +slurm` or `^mpich +slurm pmi=xxx` ( xxx is pmi or pmi2 ).
    Reference: https://slurm.schedmd.com/mpi_guide.html#mpich
-->


4. Make some packages as develop

   This step can be skiped. With this step, we can update firedrake easily in spack.

   ```bash
   spack develop py-firedrake@develop && \
   spack develop libsupermesh@develop && \
   spack develop petsc@develop && \
   spack develop slepc@develop && \
   spack develop py-fiat@develop && \
   spack develop py-finat@develop && \
   spack develop py-islpy@develop && \
   spack develop py-petsc4py@develop && \
   spack develop py-slepc4py@develop && \
   spack develop py-pyadjoint@develop && \
   spack develop py-pyop2@develop && \
   spack develop py-coffee@develop && \
   spack develop py-loopy@develop && \
   spack develop py-cgen@develop && \
   spack develop py-codepy@develop && \
   spack develop py-genpy@develop && \
   spack develop py-tsfc@develop && \
   spack develop py-ufl@develop && \
   spack develop chaco@petsc
   ```

   ````{note}
   We do not need the following package when install `int64` version:

   ```bash
   spack develop chaco@petsc
   ```
   ````

5. Concretize and install

   ```bash
   spack concretize -f 2>&1 | tee $SPACK_ENV/spack-firedrake-develop.log && \
       time spack install --fail-fast --show-log-on-error \
           --log-file $SPACK_ENV/spack-firedrake-install.log --log-format cdash
   ```

<!--
### Spack II

使用 `spack` 安装依赖包, 然后类似于 `Ubuntu` 方式安装 (需要禁用包管理器： `--no-package-manager`)

可参考 [spack-firedrake.py](./script/spack-firedrake.py).

-->

## Docker

### Images

Firedrake team provides some docker images for users to use: https://hub.docker.com/u/firedrakeproject.

For Docker on a server, there may be issues related to file permissions. Try the following script:
[`firedrake-run`](https://gist.github.com/lrtfm/93d55c6d8d5bc92d2c7551dd0048fd4f?permalink_comment_id=5024982#gistcomment-5024982)
<!--
2. lrtfm/firedrake: https://hub.docker.com/r/lrtfm/firedrake
-->

<!--
### TODO: Trimming the Docker image

The Docker image is too large, so we can consider deleting some unnecessary files.

```bash
firedrake=$HOME/firedrake
rm -rf $HOME/.cache/pip
find $firedrake -name ".git" | xargs rm -rf
find $firedrake -name "*.o" | xargs rm

rm -rf $firedrake/src/{libspatialindex,libsupermesh}

rm -rf $firedrake/src/{petsc,slepc}/src

find $firedrake -name "doc" | xargs rm -rf
find $firedrake -name "docs" | xargs rm -rf
```

### Docker commands

```bash
docker export
docker import
```
-->
