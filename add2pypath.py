#!/usr/bin/env python

import os, imp, sys

rsfpath = imp.load_source('path',
                          os.path.join(os.path.dirname(sys.argv[0]),
                                       'framework/py_pkg/path.py'))
root = os.environ.get('RSFROOT',os.getcwd())
pkgdir = os.path.split(rsfpath.get_pkgdir(root))[0]

print pkgdir
