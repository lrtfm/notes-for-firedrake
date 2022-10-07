#!/usr/bin/env python3

import logging
import subprocess
import sys
import os
import shutil
import argparse
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import shlex
import distutils.spawn

def is_tool(name):
  return distutils.spawn.find_executable(name) is not None

def log_call(cmd):
    log.info("Running command '%s'", " ".join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=os.environ)
    for line in iter(proc.stdout.readline, b''):
        log.info(line.decode().rstrip('\n'))
    ret = proc.wait()
    if ret != 0:
        e = subprocess.CalledProcessError(ret, " ".join(cmd))
        log.info(e)
        raise e

def check_call(arguments):
    try:
        log.debug("Running command '%s'", " ".join(arguments))
        log.debug(subprocess.check_output(arguments, stderr=subprocess.STDOUT, env=os.environ).decode())
    except subprocess.CalledProcessError as e:
        log.debug(e.output.decode())
        raise

def check_output(args):
    try:
        log.debug("Running command '%s'", " ".join(args))
        return subprocess.check_output(args, stderr=subprocess.STDOUT, env=os.environ).decode()
    except subprocess.CalledProcessError as e:
        log.debug(e.output.decode())
        raise

class environment(object):
    def __init__(self, **env):
        self.old = os.environ.copy()
        self.new = env

    def __enter__(self):
        os.environ.update(self.new)

    def __exit__(self, *args):
        os.environ = self.old

class workpath(object):
    def __init__(self, path):
        self.old = os.getcwd()
        self.new = path

    def __enter__(self):
        os.chdir(self.new)

    def __exit__(self, *args):
        os.chdir(self.old)

def install_spack(prefix=None):
    if prefix is None:
        prefix = os.getcwd()

    prefix = os.path.abspath(prefix)
    path = os.path.join(prefix, 'spack')

    if os.path.exists(path):
        if not os.path.isdir(path):
            log.info('Path: {} is not a directory'.format(path))
            raise Exception("Path Error", 'Path: {} is not a directory!'.format(path))
        elif os.listdir(path):
            log.info('Path: {} is not empty! Assume spack have been installed in this directory!'.format(path))
            return os.path.join(path, 'bin/spack')

    log.info("Install spack in {}".format(path))
    cmd_clone_spack = ["git", "clone", "https://github.com/spack/spack.git"]
    cmd_clone_spack.append(path)
    check_call(cmd_clone_spack)
    log.info("Spack have been installed in %s!"%prefix)

    return os.path.join(path, 'bin/spack')
    
def get_spack_env(spack, env_path):
    cmd_env_activate = [spack, "env", "activate", env_path, "--sh"]
    output = check_output(cmd_env_activate)
    output.split('\n')
    env = {}
    for line in output.split('\n'):
        if line.startswith("export"):
            pos = line.find("=")
            name, value = line[7:pos], line[pos+1:]
            env[name] = value[:-1] if value[-1] == ";" else value

    return env # spack_env

def check_create_dir(prefix):
    if not os.path.isdir(prefix):
        if os.path.exists(prefix):
            log.debug('Path %s is not a directory!'%prefix) 
            raise 
        log.info('Create path: %s'%prefix)
        os.mkdir(prefix) 
    else:
        log.debug('Path %s exists!'%prefix) 


def check_tmp_executable(spack='spack'):

    cmd = (f'{spack} python -c '
           '"import spack;'
           'import spack.util.path as spack_path;'
           'cfg = spack.config._config();'
           "bs = cfg.get_config('config')['build_stage'];"
           'tmpdir = spack_path.canonicalize_path(bs[0]);'
           'print(tmpdir, end=\'\');"')

    tempdir = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, env=os.environ).decode()

    def _check_tmp_executable():
        import tempfile
        import stat

        fid, fname = tempfile.mkstemp(dir=tempdir)
        # fid = os.open(fname, os.O_RDONLY)
        mod = os.fstat(fid).st_mode
        os.fchmod(fid, mod | stat.S_IXUSR)
        os.close(fid)
        ret = os.access(fname, os.X_OK)
        os.remove(fname)
        return ret

    ret = _check_tmp_executable()
    return ret, tempdir



# if __name__ == "__main__":

install_script_url = 'https://raw.githubusercontent.com/firedrakeproject/firedrake/%s/scripts/firedrake-install'
# 'https://raw.githubusercontent.com/firedrakeproject/firedrake/JDBetteridge/python3.10/scripts/firedrake-install'

parser = argparse.ArgumentParser(description='Install firedrake based on spack')
parser.add_argument('--logfile', type=str, help='output log file name', default='spack-firedrake-install.log')
parser.add_argument('--prefix', type=str, help='Where to install spack env. Default: $HOME/opt', default=None)
parser.add_argument('--complex', action='store_true', help='Install complex version', default=False)
parser.add_argument('--int64', action='store_true', help='Install int64 version', default=False)
parser.add_argument('--debug', action='store_true', help='Install PETSc with --with-debugging option', default=False)
parser.add_argument('--firedrake-branch', type=str, help='path or url of firedrake install script!', default="master")
parser.add_argument('--skip', action='store_true', help='Skip install spack env and firedrake! (for debug)', default=False)


args, unknown = parser.parse_known_args()

logfile = args.logfile
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)-6s %(message)s',
                    filename=logfile,
                    filemode='w')

# Log to console at INFO level
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
console.setFormatter(formatter)
logging.getLogger().addHandler(console)
log = logging.getLogger()

log.info('Options for this script: %s'%args)
log.info('Options will pass to firedrake install script: %s'%unknown)

if args.prefix is None:
    HOME = os.environ['HOME']
    prefix = "{HOME}/opt".format(HOME=HOME)
else:
    prefix = args.prefix

prefix = os.path.abspath(prefix)
log.info('The prefix is %s'%prefix)

check_create_dir(prefix)
need_install_spack = False
spack = distutils.spawn.find_executable("spack")
if spack is None:
    log.debug('Spack not founded in $PATH! Will check the default install path or install it!')
    need_install_spack = True
    spack = install_spack(prefix)

log.info('Spack founded: ' + spack)

spack_env_path = "{prefix}/firedrake-env".format(prefix=prefix)
spack_env_config = """
# add package specs to the `specs` list
# This is a Spack Environment file.
#
# It describes a set of packages to be installed, along with
# configuration settings.
spack:
# add package specs to the `specs` list
  specs: [openblas, git, autoconf, automake, libtool, cmake, zlib, bison, pkgconf,
    py-pip, boost, py-virtualenv, python+tkinter, flex@2.6.4, gmake, tmux, mesa-glu,
    scorep,gdb]
  packages:
    all:
      compiler: [gcc@9.4.0]
  view: true
  concretizer:
    unify:
      true
"""

ret, tempdir = check_tmp_executable(spack)
if ret:
    log.info(f'The files in {tempdir} is executable.')
else:
    _msg = '''
We will set `build_stage` to '$user_cache_path/stage' for the firedrake env.
You can edit the file `~/.spack/config.yaml` for spack globally:

    ```
    config:
        build_stage:
        - $user_cache_path/stage
    ```
'''
    log.info(f'''The files in {tempdir} is not executable.''')
    log.info(_msg)

    spack_env_config += """
  config:
    build_stage:
      - $user_cache_path/stage
"""

check_create_dir(spack_env_path)
spack_env_file = os.path.join(spack_env_path, 'spack.yaml')

with open(spack_env_file, 'w') as f:
    f.write(spack_env_config)

spack_env = get_spack_env(spack, spack_env_path)
log.debug("Environment variables")
log.debug(os.environ)
log.debug("Environment variables of spack env before `spack install`")
log.debug(spack_env)

with environment(**spack_env):
    if not args.skip:
        log_call([spack, "concretize", "-f"])
        log_call([spack, "install"])
    else:
        log.debug("Skip install the env in spack")

spack_env = get_spack_env(spack, spack_env_path)

install_script_url = install_script_url%args.firedrake_branch
with workpath(spack_env_path):
    check_call(['curl', '-O', install_script_url])

spack_env["PETSC_CONFIGURE_OPTIONS"] = os.environ.get("PETSC_CONFIGURE_OPTIONS", '')
install_cmd = ['python3', 'firedrake-install', '--no-package-manager', '--disable-ssh']

venv_name = "firedrake"
if args.complex:
    install_cmd += ["--complex"]
    venv_name = "firedrake-complex"

if args.int64:
    install_cmd += ["--petsc-int-type", "int64"]
    venv_name += "-int64"
    spack_env["PETSC_CONFIGURE_OPTIONS"] += '--download-scalapack --download-mumps'

if args.debug:
    spack_env["PETSC_CONFIGURE_OPTIONS"] += '--with-debugging'

install_cmd += ["--slepc",
                "--with-blas=download",
                "--pip-install", "tqdm", 
                "--pip-install", "meshio", 
                "--pip-install", "gmsh",
                # "--pip-install", "ipywidgets",
               ]

if args.firedrake_branch != "master":
    found = False
    for i, name in enumerate(unknown):
        if name == '--package-branch' and unknown[i+1] == 'firedrake':
            found = True
            break

    if not found:
        install_cmd += ['--package-branch', 'firedrake', args.firedrake_branch]

install_cmd += unknown

found = False
for i, name in enumerate(install_cmd):
    if "--venv-name" == name:
        venv_name = install_cmd[i+1]
        found = True

if not found:
        install_cmd += ['--venv-name', venv_name]

spack_env.pop("PYTHONPATH", None)
log.debug("Environment variables of spack env after `spack install`")
log.debug(spack_env)
with environment(**spack_env), workpath(spack_env_path):
    log.debug("Environment variables in block for installing firedrake")
    log.debug(os.environ)
    if not args.skip:
        log_call(install_cmd)

start_python_in_env = """
#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
FIRDRAKE_BASENAME=$( basename $( dirname $SCRIPT_DIR ) )  # firedrake
MY_SPACK_DIR=$( dirname $( dirname $SCRIPT_DIR ) )           # /../firedrake-env
MY_SPACK_LOG=$(spack env status)

if [ "`basename \"$VIRTUAL_ENV\"`" != "$FIRDRAKE_BASENAME" ]; then
        if [ "$VIRTURE_ENV" != "" ]; then
                deactivate
        fi
        if `echo "$MY_SPACK_LOG" | grep -v "$MY_SPACK_DIR" > /dev/null 2>&1`; then
                if `echo "$MY_SPACK_LOG" | grep -v "No" > /dev/null 2>&1`; then
                        spack env deactivate
                fi
                spack env activate $MY_SPACK_DIR -p
                unset PYTHONPATH
        fi
        . $SCRIPT_DIR/activate
fi

export MPIR_CVAR_ENABLE_GPU=0
export OMP_NUM_THREADS=1
# pass all the paramter to pthon
$SCRIPT_DIR/python $@

"""

start_python_script = os.path.join(spack_env_path, venv_name, 'bin', 'python_in_env')
with open(start_python_script, 'w') as f:
    f.write(start_python_in_env)

check_call(['chmod', '+x', start_python_script])

activate_fun = """
# --- FIREDRAKE FUN BEGIN --- #
function activate_firedrake() {{
  FD_SPACK_ENV_PATH=$1
  FD_ENV_NAME=$2
  {spack} env activate $FD_SPACK_ENV_PATH -p
  unset PYTHONPATH
  # export MPIR_CVAR_ENABLE_GPU=0
  export OMP_NUM_THREADS=1
  source $FD_SPACK_ENV_PATH/$FD_ENV_NAME/bin/activate
  echo "Firedrake in $FD_SPACK_ENV_PATH/$FD_ENV_NAME activated!"
}}

function activate_petsc() {{
  FD_SPACK_ENV_PATH=$1
  FD_ENV_NAME=$2
  {spack} env activate $FD_SPACK_ENV_PATH -p
  # export MPIR_CVAR_ENABLE_GPU=0
  export PETSC_DIR=$FD_SPACK_ENV_PATH/$FD_ENV_NAME/firedrake/src/petsc
  export PETSC_ARCH=default
  export PATH=\$PETSC_DIR/\$PETSC_ARCH/bin:\$PATH
}}
# --- FIREDRAKE FUN END --- #
"""

alias_template = """
# --- FIREDRAKE ALIAS BEGIN --- #
alias {venv_name}="activate_firedrake '{spack_env_path}' '{venv_name}'"
alias {venv_name}-petsc="activate_petsc '{spack_env_path}' '{venv_name}'"
# --- FIREDRAKE ALIAS END --- #
"""

alias_cmd = activate_fun.format(spack=spack) + alias_template.format(spack_env_path=spack_env_path, venv_name=venv_name)

activate_script = os.path.join(spack_env_path, venv_name, 'bin', 'firedrake-env.sh')
with open(activate_script, 'w') as f:
    f.write(alias_cmd)

log.info("")
log.info("*"*80)
if need_install_spack:
    log.info("")
    log.info("As spack is installed but not in PATH, you should add the follwing code to")
    log.info(" `.bashrc`:")
    log.info("")
    log.info("\tsource {prefix}/spack/share/spack/setup-env.sh".format(prefix=prefix)) 
    log.info("")
log.info("The following code can be added to .bashrc. Then you can activate")
log.info("firedrake by `{venv_name}` after restart the terminal".format(venv_name=venv_name))
log.info("")
log.info("\tsource " + activate_script)
log.info("")
log.info("Enjoy!")
log.info("*"*80)
log.info("")

