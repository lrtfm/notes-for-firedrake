# [Installation of Firedake](https://www.firedrakeproject.org/download.html)


<!--
> {sub-ref}`today` | {sub-ref}`wordcount-words` words | {sub-ref}`wordcount-minutes` min read

{attribution="Hamlet act 4, Scene 5"}
> We know what we are, but know not what we may be.
-->

To install firedrake, the computer should have access to the Internet.
Otherwise, please refer to {doc}`install_without_internet_access`. 

## Ubuntu

The easiest way to intall firedrake is to download the installation script `firedrake-install` and run it using Python.
This method will intall the real number version by defaults.

1. Download the installation script

   ```bash
   curl -O \
   https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
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

### Installation `real-int32` and/or `real-int32-debug`

1. Download the installation script

   ```bash
   curl -O \
   https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   ```

2. Enable PETSc's debug option (optional)

   ```bash
   DEBUG='-debug'
   sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' firedrake-install
   ```

3. Update the package of the system

   ```bash
   sudo apt-get update
   sudo apt-get install pkg-config # for p4est
   ```

4. Install

   ```bash
   PETSC_CONFIGURE_OPTIONS=" \
       --download-fftw --download-mmg \
       --download-p4est --download-parmmg --download-triangle \
       --download-tetgen --download-ctetgen --download-hpddm --download-libpng \
       --download-slepc --download-pragmatic --download-eigen" \
   python3 firedrake-install --disable-ssh \
       --documentation-dependencies \
       --venv-name $HOME/firedrake/real-int32$DEBUG
   ```

### Installation `complex-int32` and/or `complex-int32-debug`

1. Download the installation script

   ```bash
   curl -O \
   https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   ```

2. Enable PETSc's debug option (optional)

   ```bash
   DEBUG='-debug'
   sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' firedrake-install
   ```

3. Update the package of the system

   ```bash
   sudo apt-get update
   sudo apt-get install pkg-config # for p4est
   ```

4. Install

   ```bash
   PETSC_CONFIGURE_OPTIONS=" \
       --download-fftw --download-mmg \
       --download-p4est --download-parmmg --download-triangle \
       --download-tetgen --download-ctetgen --download-hpddm --download-libpng \
       --download-slepc --download-pragmatic --download-eigen \
       --download-scalapack --download-mumps" \
   python3 firedrake-install --disable-ssh \
       --documentation-dependencies  \
       --petsc-int-type int32 --complex \
       --venv-name $HOME/opt/firedrake/complex-int32$DEBUG
    ```

###  Install `complex-int64` and/or `complex-int64-debug`

1. Download the installation script

   ```bash
   curl -O \
   https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
   ```

2. Enable PETSc's debug option (optional)

   ```bash
   DEBUG='-debug'
   sed -i.bak -e 's/\(--with-debugging=\)0/\11/g' firedrake-install
   ```

3. Update the package of the system

   ```bash
   sudo apt-get update
   sudo apt-get install pkg-config
   ```

4. Install

   ```bash
   PETSC_CONFIGURE_OPTIONS=" \
       --download-fftw --download-mmg \
       --download-p4est --download-parmmg --download-triangle \
       --download-tetgen --download-ctetgen --download-hpddm --download-libpng \
       --download-slepc --download-scalapack --download-mumps" \
   python3 firedrake-install --disable-ssh \
       --documentation-dependencies  \
       --petsc-int-type int64 --complex \
       --venv-name $HOME/firedrake/complex-int64$DEBUG
   ```

   ```{note}
   `pragmatic` cannot be used with `int64`
   ```

### Installation Example with MKL

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

2. Update the packages of the system

   ```bash
   sudo apt-get update
   sudo apt-get install pkg-config # for p4est
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
       --venv-name firedrake/real-int32-mkl-debug
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

### Some notes on petsc

#### PETSc with X

1. Install `libx11-dev`

   ```bash
   sudo apt install libx11-dev
   ```

2. Add `--with-x=1` to `PETSC_CONFIGURE_OPTIONS`, and then follow the instructions in the previous section.

#### Add petsc bin to path

We can define command `add-petsc-bin`. Executing it in activated firedrake env will add the petsc/bin to PATH

```bash
alias add-petsc-bin='export \
    PATH=$PATH:$(dirname $(which python))/../src/petsc/lib/petsc/bin:$(\
    dirname $(which python))/../src/petsc/default/bin'

alias firedrake-mkl="export OMP_NUM_THREADS=1 && \
    source ~/firedrake/real-int32-mkl-debug/bin/activate && add-petsc-bin"
```

#### Download package for petsc

Sometimes, some of the packages that petsc depends on cannot be downloaded automatically.
We can add the option 

```bash
--with-packages-download-dir=<path/to/petsc/packages>
```

to obtain the list of required packages,
and then download these packages manually and put them into the path. Afterwards, configure it again with the above option.

The following python script can be used to download multiple packages.
Please modify the corresponding commands according to your needs.

```python
packages = {
# "scalapack": ['git://https://github.com/Reference-ScaLAPACK/scalapack', 
#               'https://github.com/Reference-ScaLAPACK/scalapack/archive/5bad7487f496c811192334640ce4d3fc5f88144b.tar.gz'],
"pastix": ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/pastix_5.2.3.tar.bz2'],
}
fail = {}
for name, paths in packages.items():
    print(name)
    flag = False
    for path in paths:
        print(f'try path: {path}')
        if path.startswith('git'):
            ret = os.system(f'git clone {path[6:]}')
        else:
            ret = os.system(f'curl -L -x socks5h://localhost:5000 -O {path}')
        if ret == 0:
            flag = True
            break

    if flag == False:
        fail[name] = paths
        print(f'Fail to download {name}: {paths}')

print('packages failed to download:')
print(fail)
```
### Test

```bash
source firedrake/bin/activate
cd $VIRTUAL_ENV/src/firedrake
pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
```
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

    Generate config file:
    
    ```bash
    jupyter notebook --generate-config
    ```

    Set `use_redirect_file` to `False` in file `~/.jupyter/jupyter_notebook_config.py`
    
    ```python
    c.NotebookApp.use_redirect_file = False
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

        Now you need check the python path in `kernel.json`. Make sure it is the python in your env. Otherwise, correct it.


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

可参考如下脚本:

https://raw.githubusercontent.com/lrtfm/notes-for-firedrake/main/scripts/spack-firedrake.py
-->

## Docker

### Images

1. firedrake team: https://hub.docker.com/u/firedrakeproject.

2. lrtfm/firedrake: https://hub.docker.com/r/lrtfm/firedrake

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