# Debug with tmux and tmux-mpi

1. install gdb tmux

```bash
sudo apt update && sudo apt -y install gdb tmux
```

2. install tmux-mpi

```bash
git clone https://github.com/crigler/dtach
cd dtach
./configure
make
mkdir -p $HOME/bin
cp dtach $HOME/bin
export PATH=$PATH:$HOME/bin
pip install --upgrade --no-cache-dir git+https://github.com/wrs20/tmux-mpi@master
```

3. 

```bash
tmux-mpi 16 gdb -ex run --args $(which python) $HOME/hsolver/hsolver.py -config /work/tests/model02/car_cabin_exact.yml -meshfile meshes/mesh1.msh -degree 2 -cip_enable true -freqs 1000:2000:500
```

```bash
#!/bin/bash

mpiexec -n 16 \
python3 $HOME/hsolver/hsolver.py \
        -config /work/tests/model02/car_cabin_exact.yml \
        -meshfile meshes/mesh1.msh \
        -degree 2 \
        -cip_enable true \
        -freqs 1000:2000:500
```


```bash
Thread 1 "python" received signal SIGBUS, Bus error.
__memmove_evex_unaligned_erms () at ../sysdeps/x86_64/multiarch/memmove-vec-unaligned-erms.S:708
708     ../sysdeps/x86_64/multiarch/memmove-vec-unaligned-erms.S: No such file or directory.
(gdb) bt
#0  __memmove_evex_unaligned_erms ()
    at ../sysdeps/x86_64/multiarch/memmove-vec-unaligned-erms.S:708
#1  0x00007feb09870bfb in memcpy (__len=47424, __src=<optimized out>,
    __dest=<optimized out>)
    at /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29
#2  MPID_nem_mpich_sendv_header (again=<synthetic pointer>, vc=0x55a9d0d49bb0,
    n_iov=<synthetic pointer>, iov=<synthetic pointer>)
    at ./src/mpid/ch3/channels/nemesis/include/mpid_nem_inline.h:396
#3  MPIDI_CH3_iSendv (vc=0x55a9d0d49bb0, sreq=<optimized out>,
    sreq@entry=0x55a9d3b5e430, iov=iov@entry=0x7ffe5ee5e0f0,
    n_iov=n_iov@entry=2) at src/mpid/ch3/channels/nemesis/src/ch3_isendv.c:69
#4  0x00007feb0985c7dc in MPIDI_CH3_EagerContigIsend (
    sreq_p=sreq_p@entry=0x7ffe5ee5e268,
    reqtype=reqtype@entry=MPIDI_CH3_PKT_EAGER_SEND, buf=<optimized out>,
    data_sz=<optimized out>, rank=rank@entry=2, tag=<optimized out>,
    comm=0x55a9d327aab0, context_offset=0) at src/mpid/ch3/src/ch3u_eager.c:545
    ...
(gdb) disassemble
Dump of assembler code for function __memmove_evex_unaligned_erms:
   0x00007feb0c349f40 <+0>:     endbr64
   0x00007feb0c349f44 <+4>:     mov    %rdi,%rax
   0x00007feb0c349f47 <+7>:     cmp    $0x20,%rdx
   0x00007feb0c349f4b <+11>:    jb     0x7feb0c349f80 <__memmove_evex_unaligned_erms+64>
   0x00007feb0c349f4d <+13>:    vmovdqu64 (%rsi),%ymm16
   0x00007feb0c349f53 <+19>:    cmp    $0x40,%rdx
   0x00007feb0c349f57 <+23>:    ja     0x7feb0c34a000 <__memmove_evex_unaligned_erms+192>
    ...
   0x00007feb0c34a286 <+838>:   add    %rdi,%rsi
   0x00007feb0c34a289 <+841>:   sub    %rdi,%rcx
=> 0x00007feb0c34a28c <+844>:   rep movsb %ds:(%rsi),%es:(%rdi)
   0x00007feb0c34a28e <+846>:   vmovdqu64 %ymm16,(%r8)
   0x00007feb0c34a294 <+852>:   vmovdqu64 %ymm17,0x20(%r8)
    ...
--Type <RET> for more, q to quit, c to continue without paging--q
Quit
(gdb) p/x $rsi
$13 = 0x7fe9f43c93a0
(gdb) p/x $rdi
$14 = 0x7feb016a4000
(gdb) p/x $rcx
$15 = 0x25b0
(gdb) p $rcx
$16 = 9648
(gdb) p/x $rdi + $rcx
$17 = 0x7feb016a65b0
(gdb) p/x $rsi
$13 = 0x7fe9f43c93a0
(gdb) p/x $rdi
$14 = 0x7feb016a4000
(gdb) info proc mapping
      ...
      0x7feafe633000     0x7feafe634000     0x1000     0x6000  r--p   /usr/lib/python3.10/lib-dynload/_bz2.cpython-310-x86_64-linux-gnu.so
      0x7feafe634000     0x7feafe635000     0x1000     0x7000  rw-p   /usr/lib/python3.10/lib-dynload/_bz2.cpython-310-x86_64-linux-gnu.so
      0x7feafe635000     0x7feafebf7000   0x5c2000        0x0  rw-p
      0x7feafebf7000     0x7feb03afc000  0x4f05000        0x0  rw-s   /dev/shm/mpich_shar_tmpMCilHL (deleted)
      0x7feb03afc000     0x7feb03bfc000   0x100000        0x0  rw-p
      ...
```

```bash
firedrake@985e250473ea:~$ df -h
Filesystem      Size  Used Avail Use% Mounted on
overlay         102T  467G   96T   1% /
tmpfs            64M     0   64M   0% /dev
shm              64M   64M     0 100% /dev/shm
/dev/sda1       102T  467G   96T   1% /work
tmpfs           504G     0  504G   0% /proc/asound
tmpfs           504G     0  504G   0% /proc/acpi
tmpfs           504G     0  504G   0% /proc/scsi
tmpfs           504G     0  504G   0% /sys/firmware
```


```bash
docker run --rm -it --name ubuntu --shm-size=2gb ubuntu
```

```bash
time bash -c 'for i in `seq 20`; \
    do MPICH_NO_LOCAL=1 \
    hsolver -n 16 \
    -config /work/tests/model02/car_cabin_exact.yml \
    -meshfile meshes/mesh1.msh \
    -degree 2 \
    -cip_enable true \
    -freqs 1000:2000:500; done'

Total elapsed time: 138.46 secs.

real    47m26.238s
user    496m5.754s
sys     261m26.813s
firedrake@8391e6c501ca:~$ exit
exit
```

```bash
time bash -c 'for i in `seq 20`; \
    do hsolver -n 16 \
    -config /work/tests/model02/car_cabin_exact.yml \
    -meshfile meshes/mesh1.msh \
    -degree 2 \
    -cip_enable true \
    -freqs 1000:2000:500; done'


Total elapsed time: 137.61 secs.

real    47m15.103s
user    608m18.470s
sys     146m2.341s
```

```bash
mpiexec -n 16 -bind-to core -map-by socket python /path/to/script.py
```
