#!/usr/bin/env python

import os, imp

rsfpath = imp.load_source('path','framework/py_pkg/path.py')
root = os.environ.get('RSFROOT',os.getcwd())
pkgdir = os.path.split(rsfpath.get_pkgdir(root))[0]

print pkgdir
