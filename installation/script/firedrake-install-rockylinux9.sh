#!/bin/bash
# 
# Install firedrake on Rocky Linux 9

set -e

# export http_proxy=http://127.0.0.1:1087;export https_proxy=http://127.0.0.1:1087;export ALL_PROXY=socks5://127.0.0.1:1080

sudo dnf config-manager --enable crb # ninja-build in crb
sudo dnf groupinstall --with-optional 'Development Tools'
sudo dnf install gcc-gfortran libcurl-devel libxml2-devel \
    zlib-devel boost-devel ninja-build \
    python3-devel python3-pip python3-tkinter

VENV_NAME=${VENV_NAME:-firedrake-real-int32}
ACTIVATE_FILE=$VENV_NAME/bin/activate

curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
python3 firedrake-install --no-package-manager --venv-name $VENV_NAME --with-blas=download

# Install some packages
source $ACTIVATE_FILE
pip install jupyterlab gmsh tqdm
deactivate

exit 0 # Exit

# Patch the activate file for users other than the owner
cp $ACTIVATE_FILE $ACTIVATE_FILE.$(date +"%Y%m%d%H%M%S").bak

cat <<'EOF' >> $ACTIVATE_FILE
if [ $(id -u) -ne $(id -u $(stat -c '%U' $(which python))) ] ; then
    if [ -z "${PYOP2_CACHE_DIR:-}" ]; then
        PYOP2_CACHE_DIR=$HOME/.cache/pyop2
        export PYOP2_CACHE_DIR
    fi
    
    if [ -z "${FIREDRAKE_TSFC_KERNEL_CACHE_DIR:-}" ]; then
        FIREDRAKE_TSFC_KERNEL_CACHE_DIR=$HOME/.cache/tsfc
        export FIREDRAKE_TSFC_KERNEL_CACHE_DIR
    fi
fi
EOF

cat <<'EOF' >> $ACTIVATE_FILE
if [ $(id -u) -ne $(id -u $(stat -c '%U' $(which python))) ] ; then
    if [ -n "${PIP_PREFIX:-}" ] ; then
        _OLD_PIP_PREFIX="${PIP_PREFIX:-}"
    fi

    PIP_PREFIX="$HOME/.local/share/$(basename $VIRTUAL_ENV)-$(echo $(which python) | sha256sum | awk '{print $1}')"
    export PIP_PREFIX

    if [ -n "${PYTHONPATH:-}" ] ; then
        _OLD_PYTHONPATH="${PYTHONPATH:-}"
    fi

    PYTHONPATH=$(python3 -c "
import os
import sysconfig

print(\":\".join([
    os.path.join(\"$PIP_PREFIX\", os.path.relpath(sysconfig.get_paths()[_], start=os.sys.prefix))
    for _ in (\"platlib\", \"purelib\")]))
")
    export PYTHONPATH

    PATH="$PIP_PREFIX/bin:$PATH"
    export PATH
fi
EOF

sed -i '/deactivate () {/a\
    # reset old environment variables
    if [ -n "${_OLD_PIP_PREFIX:-}" ] ; then\
        PIP_PREFIX="${_OLD_PIP_PREFIX:-}"\
        export PIP_PREFIX\
        unset _OLD_PIP_PREFIX\
    fi\
    \
    if [ -n "${_OLD_PYTHONPATH:-}" ] ; then\
        PYTHONPATH="${_OLD_PYTHONPATH:-}"\
        export PYTHONPATH\
        unset _OLD_PYTHONPATH\
    fi\
' $ACTIVATE_FILE

