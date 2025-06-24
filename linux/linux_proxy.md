# Proxy for Linux 

1. Enable proxy through dynamic port forwarding in ssh (the sockets proxy port is 5000)

    ```bash
    # Replace 'user' with your username and 'hostname' with the remote host address
    ssh -vv -N -D 5000 user@hostname
    ```

    Then you can set environment variables to route all your terminal traffic through the proxy.

    ```bash
    export http_proxy="socks5h://127.0.0.1:5000"
    export https_proxy="socks5h://127.0.0.1:5000"
    ```

    This will ensure that most command-line tools (such as `wget`, `curl`, `git`, etc.) use the SOCKS5 proxy by default.

1. apt proxy 

    Assume the proxy is on localhost at port 5000, using the SOCKS5 protocol.

    1. `-o` option

        ```bash
        sudo apt -o Acquire::http::proxy="socks5h://127.0.0.1:5000" update
        ```

    2. configuration file `/etc/apt/apt.conf`

        ```
        Acquire::http::proxy "socks5h://127.0.0.1:5000/";
        Acquire::ftp::proxy "socks5h://127.0.0.1:5000/";
        Acquire::https::proxy "socks5h://127.0.0.1:5000/";
        ```

1. curl proxy

   Use the configuration file `~/.curlrc` or the command line option `-x`

    ```bash
    curl -x socks5h://localhost:5000 -L -O <url>
    ```