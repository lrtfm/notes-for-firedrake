#!/bin/env python3
import os
import sys


if len(sys.argv) < 2:
    print(f'Usage: {sys.argv[0]} <script>')
    exit(0)

scripts = sys.argv[1]

cmd = '''env -i /bin/bash --noprofile --norc <<EOF
source {scripts} >> /dev/null
env | grep -v "^PWD=\|^SHLVL=\|^_="
EOF
'''

cmd_path = "env -i /bin/bash --norc --noprofile -c 'echo $PATH'" 
pathes = os.popen(cmd_path).read()

PATH_REMOVE = pathes.split(':')

def get_need_process_envs(content):
    envs = {}
    unset = []
    for line in content.splitlines():
        name_value = line.split('=')
        if len(name_value) == 1:
            name = name_value[0]
            value = ''
        else:
            name, value = name_value
        if ':' in value:
            envs[name] = value.split(':')
        else:
            unset.append(name)

#    PATH_REMOVE = ['/usr/bin', '/bin', '/usr/local/bin']

    for value in PATH_REMOVE:
        if value in envs['PATH']:
            envs['PATH'].remove(value)

    return unset, envs

content = os.popen(cmd.format(scripts=scripts)).read()

print(content)
print('\n'*5)

unset, envs = get_need_process_envs(content)

for name, value in envs.items():
    if name not in os.environ:
        continue
    orig = os.environ[name].split(':')
    # print(name, value, orig)
    for v in value:
        while v in orig:
            orig.remove(v)
    while '' in orig:
        orig.remove('')
    if len(orig) > 0:
        print(f"{name}={':'.join(orig)}; export {name}")
    else:
        print(f"unset {name}")

for name in unset:
    print(f'unset {name}')

