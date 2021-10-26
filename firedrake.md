# Firedrake

## Parallel computing

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