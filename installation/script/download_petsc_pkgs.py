#!/usr/bin/env python3

import os
import sys
import subprocess
from textwrap import indent


class workspace(object):
    def __init__(self, path):
        self.old = os.getcwd()
        self.new = path

    def __enter__(self):
        os.chdir(self.new)

    def __exit__(self, *args):
        os.chdir(self.old)


def read_pkgs_info(file):
    packages = {}
    for line in file:
        index = line.find('[')
        if index < 0:
            continue
        name = line[:index].strip()
        packages[name] = eval(line[index:])
    return packages


def load_pkgs_info(filename=None):
    if filename is None:
        return read_pkgs_info(sys.stdin)

    with open(filename, 'r') as file:
        return read_pkgs_info(file)


def download(pkgs_info):
    fail = {}
    n_pkgs = len(pkgs_info)
    print(f"Downloading {n_pkgs} packages to {os.getcwd()}...")
    for i, (name, paths) in enumerate(pkgs_info.items()):
        flag = False
        n_urls = len(paths)
        prefix = f"[{i+1:02d}/{n_pkgs:02d}]"
        print(f"{prefix} Download {name}...")
        for j, path in enumerate(paths):
            if path.startswith("git"):
                cmd = path
            else:
                # cmd = f"curl -L -x socks5h://localhost:5000 -O {path}"
                cmd = f"curl -L -O {path}"
            print(f"  [{j+1:d}/{n_urls:d}] {cmd}")
            ret = subprocess.run(cmd, shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True)

            if ret.returncode == 0:
                flag = True
                break
            else:
                print(indent(ret.stderr, prefix=" "*4))

        if flag == False:
            fail[name] = paths
            print(f"{prefix} Failed to download {name}.")

    if len(fail) > 0:
        n_fail = len(fail)
        print(f"\n{n_fail} package{'s' if n_fail > 1 else ''} failed to download:\n")
        for name, paths in fail.items():
            print(f"{name} {paths}")
    else:
        print('\nAll packages downloaded successfully.\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Download packages.')
    parser.add_argument('-d', '--dir', default='.',
                        help='Directory to save packages.')
    parser.add_argument('filename', nargs='?',
                        help="File containing packages info,"
                             " i.e., the output of petsc configure."
                             " If not specified, read from STDIN.")

    args = parser.parse_args()
    pkgs_info = load_pkgs_info(args.filename)

    if not os.path.exists(args.dir):
        try:
            print(f"Creating directory {args.dir}...")
            os.mkdir(args.dir)
        except:
            raise Exception(f"Failed to create directory {args.dir}")

    with workspace(args.dir):
        download(pkgs_info)