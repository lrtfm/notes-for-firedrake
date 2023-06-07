# Installation without Network

If you need to install Firedrake on some HPC without internet access, you can
use the source mirror feature of spack. A mirror is a URL that points to a
directory, either on the local filesystem or on some server, containing tarballs
for all of Spack’s packages.

We assume the local host can access the network (github, etc.). If the login
node can access the network, the operations on the local host can be executed on
the login node.  Generally, HPC uses shared storage, so there is no need to
arhive and upload the downloaded packages.

In the following, we will install spack and firedrake in directory `$HOME/opt`.

Note that the multi-line commands are connected by "`&& \`" for copy and paste.

Here is a list of web pages on spack and installing firedrake using spack:

1. spack install:

    + https://spack.readthedocs.io/en/latest/getting_started.html#installation

2. spack mirror:

    + https://spack.readthedocs.io/en/latest/bootstrapping.html#creating-a-mirror-for-air-gapped-systems
    + https://spack.readthedocs.io/en/latest/mirrors.html#mirror-environment
    + https://spack.readthedocs.io/en/latest/mirrors.html#mirror-files

3. firedrake spack:

    + https://github.com/firedrakeproject/firedrake-spack
    + https://hackmd.io/@TzVnFeL0TMCb3FaAi9qYBA/ByaRskMQ5

## Local host

The following command are run on local host with internet access.

### Create installation directory

```bash
mkdir -p $HOME/opt
```

### Clone spack

```bash
cd $HOME/opt && \
git clone -c feature.manyFiles=true https://github.com/lrtfm/spack.git && \
pushd spack && \
git checkout lrtfm/develop && \
popd && \
source $HOME/opt/spack/share/spack/setup-env.sh
```

:::{note}
Here, we clone spack from `https://github.com/lrtfm/spack.git`, a fork of spack,
and use the branch `lrtfm/develop`, which have some patchs.  You could clone
spack from the offical source `https://github.com/spack/spack.git`.
:::

:::{note}
Add the following command to `$HOME/.bashrc` to enable the shell support of spack.

```bash
source $HOME/opt/spack/share/spack/setup-env.sh
```
:::

### Create mirror for bootstrap

```bash
spack bootstrap mirror --binary-packages $HOME/opt/bootstrap
```

The output looks like:

```console
==> Adding "clingo-bootstrap@spack+python %gcc target=x86_64" and dependencies to the mirror at /home/xyz/opt/bootstrap/bootstrap_cache

==> Adding "gnupg@2.3: %gcc target=x86_64" and dependencies to the mirror at /home/xyz/opt/bootstrap/bootstrap_cache
==> Adding "patchelf@0.13.1: %gcc target=x86_64" and dependencies to the mirror at /home/xyz/opt/bootstrap/bootstrap_cache
==> Adding "gnuconfig" and dependencies to the mirror at /home/xyz/opt/bootstrap/bootstrap_cache
==> Adding binary packages from "https://github.com/spack/spack-bootstrap-mirrors/releases/download/v0.4/bootstrap-buildcache.tar.gz" to the mirro
r at /home/xyz/opt/bootstrap/bootstrap_cache

To register the mirror on the platform where it's supposed to be used, move "/home/xyz/opt/bootstrap" to its final location and run the following
command(s):

  % spack bootstrap add --trust local-sources <final-path>/metadata/sources
  % spack bootstrap add --trust local-binaries <final-path>/metadata/binaries
```

### Pack spack and bootstrap

```bash
tar -czvf spack.tar.gz spack
tar -czvf bootstrap.tar.gz bootstrap
```

### Clone firedrake-spack 

<!--
等待 `firedrake-spack` 修复
```bash
cd $HOME/opt && \
git clone https://github.com/firedrakeproject/firedrake-spack
```
-->

```bash
cd $HOME/opt && \
git clone https://github.com/lrtfm/firedrake-spack.git && \
pushd firedrake-spack && \
git checkout lrtfm/air-gapped-install && \
popd
```

:::{note}
Here, we clone firedrake-spack from
`https://github.com/lrtfm/firedrake-spack.git`, which have some patchs.
You may clone firedrake-spack from the offical source
`https://github.com/firedrakeproject/firedrake-spack.git`.
:::

### Pack firedrake-spack

```bash
tar -czvf firedrake-spack.tar.gz firedrake-spack
```

### Add repo to `spack` 

```bash
spack repo add firedrake-spack
```

:::{note}
This command can be run in an spack env which will be created below
:::

### Check the installation of spack

The command `spack info py-firedrake` should have the following output

```console
$ spack info py-firedrake
PythonPackage:   py-firedrake

Description:
    Firedrake is an automated system for the portable solution of partial
    differential equations using the finite element method (FEM)

Homepage: https://firedrakeproject.org

Preferred version:
    develop    [git] https://github.com/firedrakeproject/firedrake.git on branch master

Safe versions:
    develop    [git] https://github.com/firedrakeproject/firedrake.git on branch master

Deprecated versions:
    None

Variants:
    Name [Default]               When    Allowed values    Description
    =========================    ====    ==============    ===============================================

    64-bit-indices [off]         --      on, off           Install PETSc using 64bit indices
    build_system [python_pip]    --      python_pip        Build systems supported by the package
    complex [off]                --      on, off           Install Firedrake in complex mode
    minimal-petsc [off]          --      on, off           Build PETSc with minimal packages for Firedrake
    slepc [off]                  --      on, off           Install SLEPc and slepc4py

Build Dependencies:
    eigen            mpi            py-cython  py-h5py        py-mpi4py    py-pip        py-pyadjoint  py-scipy       py-sympy  py-vtk    slepc
    libspatialindex  petsc          py-fiat    py-islpy       py-numpy     py-pkgconfig  py-pyop2      py-setuptools  py-tsfc   py-wheel
    libsupermesh     py-cachetools  py-finat   py-matplotlib  py-petsc4py  py-progress   py-requests   py-slepc4py    py-ufl    python

Link Dependencies:
    eigen  libspatialindex  libsupermesh  mpi  petsc  python  slepc

Run Dependencies:
    eigen            petsc          py-finat       py-mpi4py    py-pip        py-pyop2         py-scipy       py-tsfc  slepc
    libspatialindex  py-cachetools  py-h5py        py-nbval     py-pkgconfig  py-pytest        py-setuptools  py-ufl
    libsupermesh     py-cython      py-islpy       py-numpy     py-progress   py-pytest-xdist  py-slepc4py    py-vtk
    mpi              py-fiat        py-matplotlib  py-petsc4py  py-pyadjoint  py-requests      py-sympy       python

```

Now, the contents of `$HOME/opt` should looks like this:

```console
$ ls -lha
total 214M
drwxrwxr-x  4 z2yang z2yang  112 Oct 30 15:17 .
drwxrwxr-x  3 z2yang z2yang   47 Oct 30 15:02 ..
drwxrwxr-x  5 z2yang z2yang  204 Oct 30 15:16 firedrake-spack
-rw-rw-r--  1 z2yang z2yang 211K Oct 30 15:17 firedrake-spack.tar.gz
drwxrwxr-x 10 z2yang z2yang 4.0K Oct 30 15:17 spack
-rw-rw-r--  1 z2yang z2yang 214M Oct 30 15:15 spack.tar.gz
```

## Remote host 

The following operations are execuated on compute nodes without internet access.

1. Installation commands should be run in compute nodes (Is this true?). In HPCs using `slurm`, you can use `srun` to start an interactive terminal:

   ```bash
   srun -p xahctest --pty --export=ALL -N 1 -n 64 --exclusive /bin/bash
   ```

   or use `salloc` first and then login to the nodes by using `ssh`:

   ```bash
   salloc -p xahctest -N 1 -n 4
   ```


2. You should add `slurm` as external package of spack on system using `slurm` to submit jobs:

   ```bash
   spack external find slurm
   ```

<!--
```bash
cd /home/z2yang/z2yang/server2
export HOME=`pwd`
export TERM=xterm
export PS1=' (server2) \w $ '
mkdir opt
cd opt
cp ../../local/opt/firedrake-spack.tar.gz .
cp ../../local/opt/spack.tar.gz .
```
-->

3. Requirements on compiler:

   1. As compilling `openblas@0.3.12` using `gcc@7.3.1` will result in error, we
      use `gcc@9.4.0`. Because Amazon Linux GCC 7.3.1 has the patch 
      [gcc-bug-87467](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87467),
      `spack` change the conflict rule for `openblas` when it is compiled by 
      `gcc@7` in [spack-pr-3443](https://github.com/spack/spack/pull/34443).
      However, GCC 7.3.1 on some hosts do not have this patch, which will result
      in error.

   2. The current version of petsc in firedrakeproject will break when using
      some compilers: 
      https://lists.mcs.anl.gov/pipermail/petsc-users/2023-April/048482.html.
      Patch has been added to branch `lrtfm/air-gapped-install` of
      `firedrake-spack`.

### Install spack

1. Create installation directory

   ```bash
   mkdir -p $HOME/opt
   ```

2. Upload files

   Upload `firedrake-spack.gz`, `spack.tar.gz`, and `bootstrap.tar.gz` to directory `$HOME/opt` on the server.

3. Unpack the files and install spack

   ```bash
   cd $HOME/opt && \
   tar -zxf spack.tar.gz && \
   tar -zxf bootstrap.tar.gz && \
   source $HOME/opt/spack/share/spack/setup-env.sh && \
   spack bootstrap add --trust local-sources $HOME/opt/bootstrap/metadata/sources && \
   spack bootstrap add --trust local-binaries $HOME/opt/bootstrap/metadata/binaries
   ```

   :::{note}
   Add the following command to the file $HOME/.bashrc to add shell support for spack.

   ```bash
   source $HOME/opt/spack/share/spack/setup-env.sh
   ```
   :::

   :::{tip}
   On some workstations, the content of the /tmp directory may not have execution permissions. You need to change the spack build directory as follows.

   ```bash
   mkdir -p $HOME/.spack && \
   cat > $HOME/.spack/config.yaml <<EOF
   config:
     build_stage:
       - \$user_cache_path/stage
   EOF
   ```
   :::

4. Install the `firedrake-spack` repo

   ```bash
   cd $HOME/opt && \
   tar -zxf firedrake-spack.tar.gz && \
   spack repo add firedrake-spack       # Can be run after the creation of the env
   ```

### Create spack env `firedrake`

1. Create spack env

   ```bash
   FIREDRAKE_ENV_NAME=firedrake-complex-int64 && \
   spack env create -d $FIREDRAKE_ENV_NAME && \
   spack env activate -p $FIREDRAKE_ENV_NAME && \
   spack -e $SPACK_ENV config add concretizer:unify:true
   ```

2. Add packages to the env

   You can add or delete some packages here. We will take the `complex+int64`
   version as an example.

   ```bash
   spack add python py-firedrake@develop%gcc +64-bit-indices+complex ^mpich ^openblas \
       ^petsc+mumps+scalapack+int64+complex+libyaml+parmmg+mmg ^llvm@12.0.1 \
       ^hypre+complex+int64+superlu-dist && \
   spack add py-pygmsh py-meshio py-tqdm py-pyyaml
   ```

<!--
    __Remark 1__: "On a cluster that uses slurm to submit jobs, you need to add slurm as an external package, and then change `^mpich` to `^mpich +slurm` or `^mpich +slurm pmi=xxx` ( xxx is pmi or pmi2 ).
    Reference: https://slurm.schedmd.com/mpi_guide.html#mpich
-->

3. Run `spack concretize`

    ```bash
    spack concretize -f 2>&1 | tee $SPACK_ENV/spack-firedrake-concretize.log
    ```

4. Check the directory `$SPACK_ENV`

   ```bash
   $ ls -la $SPACK_ENV
   total 620
   drwxrwxr-x 3 z2yang z2yang    118 Oct 30 16:01 .
   drwxrwxr-x 5 z2yang z2yang    147 Oct 30 15:33 ..
   drwxrwxr-x 4 z2yang z2yang     89 Oct 30 16:01 .spack-env
   -rw-rw-r-- 1 z2yang z2yang  54343 Oct 30 16:01 spack-firedrake-concretize.log
   -rw-rw-r-- 1 z2yang z2yang 572917 Oct 30 16:01 spack.lock
   -rw-rw-r-- 1 z2yang z2yang    457 Oct 30 16:01 spack.yaml
   ```

   We will need the `spack.lock` file to create mirror in local host.

### Create mirror on local host

The following commands run on __local host__.

We will create a firedrake env on local host by using the file `spack.lock`. And
then create mirror for the env. After that, we upload the mirror to remote host.

1. Download `spack.lock` from remote host to directory `$HOME/opt` on local host.

2. Create mirror (May take 10 mins)

   ```bash
   cd $HOME/opt && \
   spack env create -d firedrake-mirror-env spack.lock && \
   spack env activate -p ./firedrake-mirror-env && \
   time spack mirror create -a -d spack-firedrake-mirror 2>&1 | tee creat-mirror.logs
   ```

   The above command should have the following output:

   ```bash
   ==> Summary for mirror in file:///home/z2yang/z2yang/local/opt/spack-firedrake-mirror
   ==> Archive stats:
       0    already present
       244  added
       0    failed to fetch.

   real    10m56.048s
   user    1m1.559s
   sys     0m13.604s
   ```

   If there are some fails failed to fetch, you can clean the cache first and
   then create the mirror

   ```bash
   spack clean -ds && \
   time spack mirror create -a -d spack-firedrake-mirror 2>&1 | tee creat-mirror.logs
   ```

3. Pack the mirror

   ```bash
   tar -czvf spack-firedrake-mirror.tar.gz spack-firedrake-mirror
   ```

### Add mirror

The following commands run on remote host

1. Upload mirror

   Upload `spack-firedrake-mirror.tar.gz` to directory `$HOME/opt` on the remote host.

2. Unpack the mirrors

   ```bash
   cd $HOME/opt && \
   tar -xzvf spack-firedrake-mirror.tar.gz
   ```

3. Add the mirror to spack

   ```bash
   cat > $HOME/.spack/mirrors.yaml <<EOF
   mirrors:
     local_filesystem: file://$HOME/opt/spack-firedrake-mirror
   EOF
   ```

4. Check the mirror.

   The output of `spack mirror lsit` should looks like:

   ```console
   $ spack mirror list
   local_filesystem    file://<your-home-path>/opt/spack-firedrake-mirror
   spack-public        https://mirror.spack.io
   ```

### Install Firedrake

1. Run `spack develop` to avoid some errors

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
   spack develop py-ufl@develop
   ```

   :::{note}
   The `int32` version needs the following command:
   ```bash
   spack develop chaco@petsc
   ```
   :::


2. Install

   Run the following command to install (It will take 1-2 hours for the first
   time depending on the system). It may be failed. Good Luck!

   ```bash
   spack concretize -f 2>&1 | tee $SPACK_ENV/spack-firedrake-develop.log && \
   time spack install --fail-fast --show-log-on-error \
       --log-file $SPACK_ENV/spack-firedrake-install.log --log-format cdash
   ```

   The last lines of the output:

   ```console
   [+] /home/z2yang/z2yang/server2/opt/spack/opt/spack/linux-ubuntu22.04-cascadelake/gcc-11.3.0/py-firedrake-develop-il3tmyuhh37rnaww3u2yxhxcqawp3hh6
   ==> Updating view at /home/z2yang/z2yang/server2/opt/firedrake-complex-int64/.spack-env/view

   real    184m44.242s
   user    836m56.348s
   sys     97m42.901s
   ```

3. Deactivate the env

   ```bash
   despacktivate
   ```

   Feel free to ignore the following warnings

   ```console
   $ despacktivate
   ==> Warning: Skipping reversal of unreversable operation<class 'spack.util.environment.UnsetEnv'> PETSC_ARCH
   ==> Warning: Skipping reversal of unreversable operation<class 'spack.util.environment.UnsetEnv'> PETSC_ARCH
   ==> Warning: Skipping reversal of unreversable operation<class 'spack.util.environment.UnsetEnv'> PETSC_ARCH
   ==> Warning: Skipping reversal of unreversable operation<class 'spack.util.environment.UnsetEnv'> PETSC_ARCH
   ```

### Usage

1. Activate the env

   ```bash
   cd $HOME/opt && \
   spack env activate -p $FIREDRAKE_ENV_NAME
   ```

2. Test

   ```bash
   cd $SPACK_ENV/py-firedrake && \
   pytest tests/regression/ -k "poisson_strong or stokes_mini or dg_advection"
   ```