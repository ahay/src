#!/bin/env python
# Prepare a list of programs without examples for removing them from the stable version

import rsf.doc
import rsf.prog

progs = rsf.doc.progs
for prog in progs.keys():
    sfprog = rsf.doc.progs[prog]
    if not sfprog.uses:
        print sfprog.file

