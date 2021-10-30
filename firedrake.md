# Firedrake

## Installation

The doc of Firedrake [Obtaining Firedrake](https://www.firedrakeproject.org/download.html) present a detail introduction on how to install it.

### Based on spack

If you do not have access to the system's package-manager, you can install firedrak based on [`spack`](https://spack.io/).
Please see the istallation script [firedrake-space.sh](https://github.com/lrtfm/notes-for-firedrake/blob/main/firedrake-spack.sh).

## Parallel computing

It's easy to run a simulation in parallel in firedrake:
```
mpiexec -n 16 python simulation.py
```
See https://www.firedrakeproject.org/parallelism.html for more details.

<!--
### Parallel of the for loop
serial
```
items = list()
for item in items:
    process(item)
```

parallel
```
from multiprocessing.dummy import Pool as ThreadPool
items = list()
pool = ThreadPool()
pool.map(process, items)
pool.close()
pool.join()
```
-->
