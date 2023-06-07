# Proxy for Linux 
1. Enable proxy through dynamic port forwarding in ssh (the sockets proxy port is 5000)

    ```bash
    ssh -vv -N -D 5000 user@hostname
    ```

2. apt proxy 

    ```bash
    sudo apt -o Acquire::http::proxy="socks5h://127.0.0.1:5000" update
    ```

3. curl proxy

    ```bash
    curl -x socks5h://localhost:5000 -O https://url/to/you/file
    ```