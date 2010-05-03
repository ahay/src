#!/usr/bin/env python

import os, sys

if sys.version_info[:2] < (2, 7):
    import distutils.sysconfig as sysconfig
else:
    import sysconfig

root = os.environ.get('RSFROOT',os.getcwd())
std_pkgdir = sysconfig.get_python_lib()
pkgdir = std_pkgdir.replace(sysconfig.PREFIX,root,1)

print pkgdir
