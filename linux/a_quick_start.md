# Linux 快速入门

系统的学习请看 [鸟哥的私房菜](http://cn.linux.vbird.org/linux_basic/linux_basic.php),
入门推荐阅读章节请参考 [a_short_introduction_of_linux.pdf](./a_short_introduction_of_linux.pdf) 的最后几页.

## 文件系统

目录结构

```bash
/-/
  |-bin/
  |-boot/
  |-dev/
  |-etc/
  |-home/
  |     |-user1/
  |     |-user2/
  |-lib/
  |-media/
  |-mnt/
  |-opt/
  |-proc/
  |-root/
  |-run/
  |-sbin/
  |-srv/
  |-sys/
  |-tmp/
  |-usr/
  |-var/
```

## 获取帮助

通过 man 查看命令的详细信息
```
man <command>
```

man 手册有多个章节, 通过 `man <section> <command>` 来查看不同章节的命令的帮助信息.

```
1   Executable programs or shell commands
2   System calls (functions provided by the kernel)
3   Library calls (functions within program libraries)
4   Special files (usually found in /dev)
5   File formats and conventions, e.g. /etc/passwd
6   Games
7   Miscellaneous (including macro packages and conventions), e.g. man(7), groff(7)
8   System administration commands (usually only for root)
9   Kernel routines [Non standard]
```

## Shell

Shell 是用户和操作系统内核之间的交互接口, 常见的 Shell 有 Bash, zsh, csh, ksh 等. 不同的 shell 有不同的特性, 但是基本的命令是相同的.

### 环境变量 

1. PATH
2. HOME
3. PWD

4. USER
5. SHELL
6. TERM
6. LANG

### 常用命令

基本操作
+ echo
+ ls, pwd, cd, pushd, popd
+ cp, mv, rm, mkdir, touch

获取帮助
+ man, info

查看文件
+ cat, more, less, head, tail
+ tee

编辑文件
+ nano, vim, emacs

文本处理
+ cut, sort, uniq
+ grep, sed, awk
+ diff, patch

查找文件与命令
+ find, which, locate

参数
+ xargs

查看当前登录用户
+ who, whoami

shell 设置与环境变量
+ set, export, unset

修改文件权限
+ chmod, chown

查看硬盘使用情况
+ df, du

压缩与同步
+ tar
+ rsync
+ scp, sftp

下载
+ wget, curl

查看进程
+ top, htop, free 
+ ps, kill

### 重定向与管道

标准输入 (0), 标准输出 (1), 标准错误 (2)

1. `>` 重定向标准输出, `>>` 追加标准输出

    ```bash
    echo 'abc' > a.txt
    cat a.txt
    echo 'abc' >> a.txt
    cat a.txt
    ```

2. `2>` 重定向标准错误, `2>>` 追加标准错误

    ```bash
    cat not_a_file 2> b.txt
    cat b.txt
    cat not_b_file 2>> b.txt
    cat b.txt
    ```

    ```bash
    cat a.txt not_a_file >output.txt 2> error.txt
    ```

    ```bash
    cat a.txt not_a_file >output.txt 2>&1
    ```

3. 管道 `|`
   
    ```bash
    ls -l | grep 'abc'
    ```

4. 了解文件描述符

### 单行命令

显示内存使用情况

```bash
while true; do free -h; sleep 1; echo -en '\e[3A'; done
```

## 其他常用工具

### 终端复用
+ screen, tmux

### 版本管理
+ git

## 远程登录 Linux 服务器

### ssh

```bash
ssh <user-name>@<hostname> [-p <port>]
```

```bash
ssh -J xxx@jump-hostname yyy@hostname -p <port> 
```

### ssh 配置文件

通过中间服务器登录

```bash
cat >> ~/.ssh/config << EOF
Host workstation
    Hostname ip-address-or-hostname
    Port 22
    User username
    IdentityFile ~/.ssh/id_ed25519
    ProxyJump someone@proxy.com
EOF
```

```bash
ssh workstation
```

### ssh 免密登录设置

1. 生成密钥对

    ```bash
    ssh-keygen -t ed25519 -C "your-customer-comment"
    ```

    ```bash
    ssh-keygen -t rsa -b 4096 -C "your-customer-comment"
    ```


2. 拷贝公钥到服务器

    ```bash
    ssh-copy-id <user-name>@<hostname>
    ```

### Windows 客户端:

1. 客户端 MobaXterm  https://mobaxterm.mobatek.net

2. Putty

## `grep` 和正则表达式

一下内容摘自 grep 手册

```text
REGULAR EXPRESSIONS
       A  regular  expression  is  a  pattern  that  describes a set of strings.
       Regular   expressions   are   constructed   analogously   to   arithmetic
       expressions, by using various operators to combine smaller expressions.

       grep  understands  three different versions of regular expression syntax:
       “basic” (BRE), “extended” (ERE) and “perl” (PCRE).  In GNU grep there  is
       no  difference  in  available  functionality  between  basic and extended
       syntaxes.  In other implementations, basic regular expressions  are  less
       powerful.    The   following  description  applies  to  extended  regular
       expressions; differences for basic  regular  expressions  are  summarized
       afterwards.    Perl-compatible   regular   expressions   give  additional
       functionality, and are documented in  pcresyntax(3)  and  pcrepattern(3),
       but work only if PCRE is available in the system.

       The  fundamental building blocks are the regular expressions that match a
       single character.  Most characters, including all letters and digits, are
       regular  expressions  that  match  themselves.   Any  meta-character with
       special meaning may be quoted by preceding it with a backslash.

       The period . matches any single character.  It is unspecified whether  it
       matches an encoding error.

   Character Classes and Bracket Expressions
       A  bracket  expression  is  a list of characters enclosed by [ and ].  It
       matches any single character in that list.  If the first character of the
       list  is the caret ^ then it matches any character not in the list; it is
       unspecified whether it matches  an  encoding  error.   For  example,  the
       regular expression [0123456789] matches any single digit.

       Within   a  bracket  expression,  a  range  expression  consists  of  two
       characters separated by a hyphen.  It matches any single  character  that
       sorts between the two characters, inclusive, using the locale's collating
       sequence and character set.  For example, in the default C locale,  [a-d]
       is  equivalent  to  [abcd].   Many  locales sort characters in dictionary
       order, and in these locales [a-d] is typically not equivalent to  [abcd];
       it  might  be  equivalent  to  [aBbCcDd],  for  example.   To  obtain the
       traditional interpretation of bracket expressions,  you  can  use  the  C
       locale by setting the LC_ALL environment variable to the value C.

       Finally,  certain  named  classes  of  characters  are  predefined within
       bracket expressions, as follows.  Their names are self  explanatory,  and
       they   are   [:alnum:],   [:alpha:],   [:blank:],  [:cntrl:],  [:digit:],
       [:graph:], [:lower:], [:print:],  [:punct:],  [:space:],  [:upper:],  and
       [:xdigit:].   For  example,  [[:alnum:]]  means  the  character  class of
       numbers and letters in the current locale.  In the  C  locale  and  ASCII
       character  set encoding, this is the same as [0-9A-Za-z].  (Note that the
       brackets in these class names are part of the symbolic names, and must be
       included  in addition to the brackets delimiting the bracket expression.)
       Most  meta-characters  lose  their   special   meaning   inside   bracket
       expressions.   To  include  a  literal  ]  place  it  first  in the list.
       Similarly, to include a literal ^ place it anywhere but first.   Finally,
       to include a literal - place it last.

   Anchoring
       The  caret  ^ and the dollar sign $ are meta-characters that respectively
       match the empty string at the beginning and end of a line.

   The Backslash Character and Special Expressions
       The symbols \<  and  \>  respectively  match  the  empty  string  at  the
       beginning  and  end of a word.  The symbol \b matches the empty string at
       the edge of a word, and \B matches the empty string provided it's not  at
       the  edge  of a word.  The symbol \w is a synonym for [_[:alnum:]] and \W
       is a synonym for [^_[:alnum:]].

   Repetition
       A regular expression  may  be  followed  by  one  of  several  repetition
       operators:
       ?      The preceding item is optional and matched at most once.
       *      The preceding item will be matched zero or more times.
       +      The preceding item will be matched one or more times.
       {n}    The preceding item is matched exactly n times.
       {n,}   The preceding item is matched n or more times.
       {,m}   The  preceding  item  is  matched  at most m times.  This is a GNU
              extension.
       {n,m}  The preceding item is matched at least n times, but not more  than
              m times.

   Concatenation
       Two  regular  expressions  may  be  concatenated;  the  resulting regular
       expression matches any string formed by concatenating two substrings that
       respectively match the concatenated expressions.

   Alternation
       Two  regular  expressions  may  be  joined  by  the infix operator |; the
       resulting regular expression matches any string matching either alternate
       expression.

   Precedence
       Repetition  takes  precedence  over  concatenation,  which  in turn takes
       precedence over alternation.  A  whole  expression  may  be  enclosed  in
       parentheses to override these precedence rules and form a subexpression.

   Back-references and Subexpressions
       The  back-reference  \n, where n is a single digit, matches the substring
       previously matched by the nth parenthesized subexpression of the  regular
       expression.

   Basic vs Extended Regular Expressions
       In  basic  regular  expressions  the meta-characters ?, +, {, |, (, and )
       lose their special meaning; instead use the backslashed versions \?,  \+,
       \{, \|, \(, and \).
```

### Basic (`-G` default) vs Extended (`-E`) Regular Expressions

In basic regular expressions the meta-characters `?`, `+`, `{`, `|`, `(`, and `)` lose their special meaning; instead use the backslashed versions `\?`,  `\+`, `\{`, `\|`, `\(`, and `\)`.

### Braket

1. Most meta-characters lose their special meaning inside bracket expressions.  
2. To include a literal `]` place it first in the list. Similarly, to include a literal `^` place it anywhere  but first.  
3. Finally, to include a literal `-` place it last.
