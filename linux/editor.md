# Editor

## VIM + vimtex in (WSL)

### 反向搜索 (从 pdf 跳转到源文件)

使用 zathura 只需简单配置：
```
let g:vimtex_view_method = 'zathura'
```

不过需要注意以下事项:
1. vim 的 *clientserver* 选项需要开启才能使用反向搜索 (backward search or inverse search).
2. vim 的 server 依赖于 X server, 所以 X server 也需要开启, vim 的server 记录在 X server 的属性中
3. X server 的属性可以通过 xprop 命令查看