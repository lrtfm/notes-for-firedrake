# Debug with tmux and tmux-mpi

1. Install tmux

    On Ubuntu:
    ```bash
    sudo apt update && sudo apt -y install tmux
    ```

    On macOS:
    ```bash
    brew install tmux
    ```

2. Install tmux-mpi

    First, install `dtach`:

    ```bash
    git clone https://github.com/crigler/dtach
    cd dtach
    ./configure
    make
    mkdir -p $HOME/bin
    cp dtach $HOME/bin
    export PATH="$HOME/bin:$PATH"
    ```

    Then, install `tmux-mpi`:

    ```bash
    pip install --upgrade --no-cache-dir git+https://github.com/wrs20/tmux-mpi@master
    ```

3. Debug with tmux-mpi

    You can use `gdb` as the debugger on Linux and Intel-based macOS systems. On macOS with Apple Silicon (arm64), use `lldb` instead of `gdb`.

    For configuration details and advanced usage, refer to the [tmux-mpi documentation](https://github.com/wrs20/tmux-mpi?tab=readme-ov-file#configuration).

    Example usage with 16 processes:

    ```bash
    tmux-mpi 16 gdb -ex run --args $(which python) test_example.py
    ```

    Replace `test_example.py` with your script as needed.


    ```console
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
