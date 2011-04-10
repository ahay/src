#!/usr/bin/env python
'''
Finds in how many projects each Madagascar program is used.
Can be run anywhere, with no arguments.
I.e. "./document.py >results.asc"
will produce something like
"sfbyte is used in 568 projects
    clip is not documented
    scalebar is not documented
    gpow is not documented
etc..."
Note: private projects are also included
'''

from rsf.doc import progs
import rsf.prog as rsfprog

def nuses(p):
    'how many times a program is used'
    n=0
    uses = progs[p].uses
    for book in uses.keys():
        for chapter in uses[book].keys():
            n = n + len(uses[book][chapter])
    return n

programs = progs.keys()
programs.sort(lambda x,y: nuses(y)-nuses(x))

for prog in programs:
    print '%s (%s) is used in %d projects' % (prog,progs[prog].file,nuses(prog))
    pars = progs[prog].pars
    for par in pars.keys():
        if not pars[par].desc.strip():
            print '\t%s is not documented' % par 
