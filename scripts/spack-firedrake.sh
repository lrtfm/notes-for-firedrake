#!/bin/bash
# A script to install firedrake and firedrake-complex on system
# where you have no access to the package-manager.
#
# We install the required packages by using spack.
# TODO: make it a python script

git clone https://github.com/spack/spack.git


echo "source `pwd`/spack/share/spack/setup-env.sh" >> ~/.bashrc

source `pwd`/spack/share/spack/setup-env.sh

FD_ENV_PATH=~/opt/firedrake-env
mkdir -p $FD_ENV_PATH
cd $FD_ENV_PATH

cat <<"EOF" > spack.yaml
# add package specs to the `specs` list
# This is a Spack Environment file.
#
# It describes a set of packages to be installed, along with
# configuration settings.
spack:
  # add package specs to the `specs` list
  specs: [openblas, git, autoconf, automake, libtool, cmake, zlib, bison, pkgconf,
    py-pip, boost, py-virtualenv, python+tkinter, flex@2.6.4, gmake, tmux, mesa-glu,
    scorep]
  packages:
    all:
      compiler: [gcc@8.3.1]
  view: true
  concretization: together
EOF

spack env activate .
echo "Install packages by spack"
spack install

curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
echo "Installing firedrake"
python3 firedrake-install --no-package-manager --disable-ssh --pip-install scipy --pip-install tqdm --venv-name firedrake
echo "Installing firedrake-complex"
python3 firedrake-install --no-package-manager --complex --disable-ssh --pip-install scipy --pip-install tqdm --venv-name firedrake-complex

echo "Writing config into .bashrc"
cat <<EOF >> ~/.bashrc
function activate_firedrake() {
        spack env activate `pwd`
        unset PYTHONPATH
        source `pwd`/firedrake/bin/activate
}
function activate_firedrake_complex() {
        spack env activate `pwd`
        unset PYTHONPATH
        source `pwd`/firedrake/bin/activate
}
alias firedrake=activate_firedrake
alias firedrake-complex=activate_firedrake_complex
EOF

echo "Done!"

