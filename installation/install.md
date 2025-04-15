# Installation of Firedake


<!--
> {sub-ref}`today` | {sub-ref}`wordcount-words` words | {sub-ref}`wordcount-minutes` min read

{attribution="Hamlet act 4, Scene 5"}
> We know what we are, but know not what we may be.
-->

See the official document [Install](https://www.firedrakeproject.org/install.html#)

Using [Docker](https://www.firedrakeproject.org/install.html#docker) with the image provided by the Firedrake project is recommended.

To install Firedrake, the computer should have access to the internet (including `github.com` and `pypi`).


## Some Notes on Installation

### Install petsc with external packages pre-downloaded into a directory

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


### Using firedrake in Jupyter-lab (which installed not in firedrake env)

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