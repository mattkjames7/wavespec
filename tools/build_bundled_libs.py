#!/usr/bin/env python3
"""Build bundled C libraries (CMake) and copy resulting .so into package data.

Run this before building wheels if you want the shared libraries included in
the packaged distribution.
"""
import os
import subprocess
import shutil
import sys

ROOT = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
DATA = os.path.join(ROOT, 'wavespec', '__data')
BUILD = os.path.join(DATA, 'build')

def run(cmd):
    print('Running:', ' '.join(cmd))
    subprocess.check_call(cmd)

def main():
    if not os.path.isdir(DATA):
        print('Could not find', DATA)
        sys.exit(1)
    try:
        run(['cmake', '-S', DATA, '-B', BUILD])
        run(['cmake', '--build', BUILD, '--config', 'Release'])
    except subprocess.CalledProcessError as e:
        print('CMake build failed:', e)
        sys.exit(2)

    candidates = [
        os.path.join(BUILD, 'liblombscargle', 'liblombscargle.so'),
        os.path.join(BUILD, 'libfilter', 'libfilter.so'),
        os.path.join(BUILD, 'liblombscargle.so'),
        os.path.join(BUILD, 'libfilter.so'),
    ]
    for src in candidates:
        if os.path.isfile(src):
            dest = os.path.join(DATA, os.path.basename(src))
            shutil.copy2(src, dest)
            os.chmod(dest, 0o755)
            print('Copied', src, '->', dest)

if __name__ == '__main__':
    main()
