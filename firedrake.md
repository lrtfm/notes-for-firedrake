# Firedrake

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
