#!/usr/bin/env python

import os, imp

rsfpath = imp.load_source('path','framework/py_pkg/path.py')
root = os.environ.get('RSFROOT',os.getcwd())
pkgdir = rsfpath.get_pkgdir(root)

print pkgdir
