# Proxy for Linux 

1. Enable proxy through dynamic port forwarding in ssh (the sockets proxy port is 5000)

    ```bash
    ssh -vv -N -D 5000 user@hostname
    ```

2. apt proxy 

    1. `-o` option

        ```bash
        sudo apt -o Acquire::http::proxy="socks5h://127.0.0.1:5000" update
        ```

    2. configure file `/etc/apt/apt.conf`

        ```
        Acquire::http::proxy "http://127.0.0.1:8000/";
        Acquire::ftp::proxy "ftp://127.0.0.1:8000/";
        Acquire::https::proxy "https://127.0.0.1:8000/";
        ```

3. curl proxy

   use configure file `~/.curlrc` or command line option `-x`

    ```bash
    curl -x socks5h://localhost:5000 -L -O <url>
    ```