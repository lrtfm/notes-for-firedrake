# Install firedrake on CentOS 7
These notes describe how to install Firedrake on a small HPC system running CentOS 7 without root privileges.
All dependencies, including GCC and Python, are installed using Spack.

> **Note:**
> 1. The login node has internet access, but direct access to `github.com` may be restricted.  
> 2. If you need to access `github.com`, use an SSH proxy as described in the [SSH proxy](#ssh-proxy) section.

## Install spack

We install spack on dir `$HOME/opt`.

Clone the repo and load spack env
```
git clone --depth=2 --branch=releases/v0.23 https://github.com/spack/spack.git $HOME/opt/spack
```

```
. $HOME/opt/spack/share/spack/setup-env.sh
```

<!--
2. Install module tool 
reference: https://spack-tutorial.readthedocs.io/en/latest/tutorial_modules.html#build-a-module-tool

```
spack install lmod
. $(spack location -i lmod)/lmod/lmod/init/bash
. $HOME/opt/spack/share/spack/setup-env.sh
```
-->

## Install dependencies

1. Install required packages:
    ```bash
    spack install gcc@9.3.0
    spack install bison cmake flex git gmake openblas ninja python tmux cuda@11.2.0 py-venv py-pip py-pysocks
    spack install openmpi +cuda cuda_arch=80,86 fabrics=ucx ^cuda@11.2.0 \
        ^ucx+cuda+verbs+cma+dc+dm+gdrcopy+mlx5_dv+thread_multiple cuda_arch=80,86
    ```

2. Create a Spack view to collect the installed tools:
    ```bash
    spack view -d no add -i $HOME/opt/tools-view \
        gcc@9.3.0 openmpi hwloc bison cmake flex git gmake ninja python tmux openblas zlib-ng ucx \
        py-venv py-pip py-pysocks
    ```

   If you want to remove a package from the view (for example, `cmake`), run:
    ```bash
    spack view -d no rm $HOME/opt/tools-view cmake
    ```

3. Add the `bin` path of view to `PATH`
    ```bash
    export PATH="$HOME/opt/tools-view/bin:$PATH"
    ```

### Notes on CUDA and Kokkos

1. **CUDA:**  
    Refer to the [CUDA 11.2 installation guide](https://docs.nvidia.com/cuda/archive/11.2.0/cuda-installation-guide-linux/index.html).  
    The GPU nodes use Nvidia driver version `460.27.04`, which is compatible with CUDA 11.2. According to the CUDA documentation, GCC versions below 10 are required.

2. **Kokkos:**  
    See the [Kokkos requirements](https://kokkos.org/kokkos-core-wiki/get-started/requirements.html).  
    Kokkos requires a compiler version greater than 8.2.0.

Based on these requirements, GCC 9.3.0 is selected. Note that GCC 9.4.0 and 9.5.0 may cause the error:  
`error: identifier "__builtin_ia32_rndscalesd_round" is undefined`


## Install Firedrake

Follow the official [Firedrake installation instructions](https://www.firedrakeproject.org/install.html).

See [install-firedrake-centos7.sh](script/install-firedake-centos7.sh) for an example.

### Configure PETSc

When building PETSc, you need to specify the locations of dependencies so PETSc can find them. Set the environment variable `$DEPS_VIEW` to your Spack view directory (e.g., `$HOME/opt/tools-view`). Add the following options to your PETSc `configure` command:

```bash
--with-hwloc-dir=$DEPS_VIEW \
--with-openblas-dir=$DEPS_VIEW \
--with-zlib-dir=$DEPS_VIEW \
```

For CUDA support, set `$CUDA_DIR` to the CUDA installation path (find it with `spack location -i cuda`) and add:

```bash
--with-cuda-dir=$CUDA_DIR
```

To enable Kokkos support, add the following options to your PETSc `configure` command (replace `80` with the compute capability of your GPU if different):

```bash
--with-cuda-arch=80 `#maybe not needed https://petsc.org/release/changes/316/#changes-3-16` \
--download-kokkos \
--download-kokkos-kernels
```

### Notes

- Since OpenMPI is built with CUDA support, you can disable CUDA at runtime by setting the environment variable:
    ```bash
    export OMPI_MCA_opal_cuda_support=0
    ```

- If running on a system without CUDA, add the CUDA stubs directory to `LD_LIBRARY_PATH` to suppress warnings about missing `libcuda.so`:
    ```bash
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$CUDA_DIR/lib64/stubs"
    ```

- If you encounter an error about `libcuda.so.1` not being found, create a symbolic link in the stubs directory:
    ```bash
    ln -s libcuda.so "$CUDA_DIR/lib64/stubs/libcuda.so.1"
    ```

## SSH proxy

This section describes how to set up an SSH SOCKS5 proxy on port 5000 and ensure it is automatically stopped when your session ends.

First, configure passwordless SSH access to your proxy server by adding the following to your `~/.ssh/config`:

```ssh
Host proxy-server
    HostName your.proxy.server.address
    User your-username
```

To start the proxy and set the necessary environment variables for `git`, `curl`, and `pip`, use the following script:

```bash
#!/bin/bash

set -e

ssh -vv -n -N -D 5000 -o ExitOnForwardFailure=yes proxy-server >ssh-proxy.log 2>&1 &
SSH_PID=$!

for i in {1..30}; do
    if ! ps -p $SSH_PID > /dev/null; then
        echo "SSH failed to start. See ssh-proxy.log for details."
        exit 1
    fi
    sleep 0.1
done

echo "Started SSH tunnel with PID $SSH_PID"

trap "echo 'Killing SSH tunnel...'; kill $SSH_PID" EXIT

export http_proxy="socks5h://localhost:5000"
export https_proxy="socks5h://localhost:5000"

# Place your commands that require proxy access below
```