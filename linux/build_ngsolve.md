# install ngsolve on ws2


1. install libffi https://www.linuxfromscratch.org/lfs/view/9.1/chapter06/libffi.html

   ```bash
   cd $HOME/pkgs
   wget https://github.com/libffi/libffi/releases/download/v3.4.4/libffi-3.4.4.tar.gz
   tar -zxvf libffi-3.4.4.tar.gz
   cd libffi-3.4.4
   ./configure --prefix=$HOME/opt/libffi-3.4.4 --disable-static --with-gcc-arch=native
   make
   make check
   make install
   ```

2. install Binutils if needed

   Ref: https://www.linuxfromscratch.org/lfs/view/9.1-systemd/chapter06/binutils.html

   ```
   tar -zxvf binutils-2.38.tar.gz
   cd binutils-2.38
   # The Binutils documentation recommends building Binutils in a dedicated build directory:
   mkdir -v build
   cd       build
   ../configure --prefix=$HOME/opt/local
   make
   make install
   ```

2. install openssl if needed

   Ref: https://www.linuxfromscratch.org/lfs/view/9.1-systemd/chapter06/openssl.html
    
    ```
    ./config --prefix=$HOME/opt/local         \
             --libdir=lib          \
             shared                \
             zlib-dynamic
    make 
    make install
    ```

3. install sqlite (optional)

    ```bash
    ./configure --prefix=$HOME/opt/local     \
                --disable-static  \
                --enable-fts5     \
                CPPFLAGS="-DSQLITE_ENABLE_FTS3=1            \
                          -DSQLITE_ENABLE_FTS4=1            \
                          -DSQLITE_ENABLE_COLUMN_METADATA=1 \
                          -DSQLITE_ENABLE_UNLOCK_NOTIFY=1   \
                          -DSQLITE_ENABLE_DBSTAT_VTAB=1     \
                          -DSQLITE_SECURE_DELETE=1          \
                          -DSQLITE_ENABLE_FTS3_TOKENIZER=1" &&
    make
    make install
    ```

4. install python

   ```bash
   cd $HOME/pkgs
   wget https://www.python.org/ftp/python/3.11.7/Python-3.11.7.tgz
   tar -zxvf Python-3.11.7.tgz
   cd Python-3.11.7
   LDFLAGS="-L$HOME/opt/local/lib -L$HOME/opt/local/lib64 -L$HOME/opt/libffi-3.4.4/lib64 \
            -Wl,-rpath=$HOME/opt/local/lib -Wl,-rpath=$HOME/opt/local/lib64 \
            -Wl,-rpath,$HOME/opt/libffi-3.4.4/lib64 -Wl,-rpath,$HOME/opt/python-3.11/lib" \
       CPPFLAGS="-I$HOME/opt/libffi-3.4.4/include -I$HOME/opt/local/include" \
       CFLAGS="-I$HOME/opt/libffi-3.4.4/include -I$HOME/opt/local/include" \
       ./configure --prefix=$HOME/opt/python-3.11 --enable-optimizations --enable-shared
   make
   make install
   ```

   Some output of make (We just ignore it. In addition, nis is deprecated):
   ```console
   The necessary bits to build these optional modules were not found:
   _curses               _curses_panel         _dbm
   _gdbm                 nis                   readline
   To find the necessary bits, look in setup.py in detect_modules() for the module's name.
   
   
   The following modules found by detect_modules() in setup.py have not
   been built, they are *disabled* by configure:
   _sqlite3
   ```

5. Add python to PATH (in .bashrc)

   ```bash
   export PATH=$HOME/opt/python-3.11/bin:$PATH
   ```

6. Install cmake

   ```bash
   wget https://github.com/Kitware/CMake/releases/download/v3.28.1/cmake-3.28.1.tar.gz
   tar -zxvf cmake-3.28.1.tar.gz
   cd cmake-3.28.1
   ./configure --prefix=$HOME/opt/cmake-3.28.1 && make -j32 && make install
   ```

7. install pybind11-stubgen

   ```bash
   python3 -m pip install pybind11-stubgen
   ```

8. install ngsolve

   ```bash
   mkdir $HOME/pkgs/ngsuite
   cd $HOME/pkgs/ngsuite
   git clone --recurse-submodules https://github.com/NGSolve/ngsolve.git ngsolve-src
   mkdir ngsolve-build
   cd ngsolve-build
   cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/ngsolve-230104 -DCMAKE_PREFIX_PATH=$HOME/opt/ngsolve-230104 -DUSE_OCC=OFF -DUSE_GUI=OFF ../ngsolve-src
   ```

   Run the following command before running `make -j32` and `make install` in the build directory.

   ```bash
   NETGENDIR=/home/21041416r/opt/ngsolve-230104/bin
   PYTHONPATH=.:/home/21041416r/opt/ngsolve-230104/lib/python3.11/site-packages
   ```

   This is to fix the error when make
   ```
   ModuleNotFoundError: No module named 'netgen'
   ```

9.  Add the following to .bashrc, which is a function to activate ngsolve

   ```bash
   activate_ngsolve () {
       export PS1="( ngsolve ) $PS1"
       export NETGENDIR=/home/21041416r/opt/ngsolve-230104/bin
       export PATH="$NETGENDIR:$PATH"
       export PYTHONPATH=".:/home/21041416r/opt/ngsolve-230104/lib/python3.11/site-packages:$PYTHONPATH"
   }
   ```

10. activate ngsolve to use

   ```bash
   activate_ngsolve
   ```

   Then you can run your ngsolve scripts

   ```
   python3 your-python-script.py
   ```