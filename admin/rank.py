#!/usr/bin/env python
'''
Establish comparative ranking of different programs and projects
'''

from rsf.doc import progs
import rsf.prog

import networkx as nx
from functools import cmp_to_key

def cmp(a, b):
    return (a > b) - (a < b) 

def output_rank(G,filename):
    rank = nx.pagerank(G)
    keys = list(rank.keys())
    keys.sort(key = cmp_to_key(lambda x,y: cmp(rank[y],rank[x])))

    ofile = open(filename,'w')
    for key in keys:
        ofile.write('%s %s\n' % (rank[key],key))
    ofile.close()
    

Gprog = nx.Graph()
Gproj = nx.Graph()

programs = list(progs.keys())

proj = {}
for p in programs:
    dirs = []
    uses = progs[p].uses
    for book in list(uses.keys()):
        for chapter in list(uses[book].keys()):
            for project in uses[book][chapter]:
                dir = '/'.join([book,chapter,project])
                dirs.append(dir)

                if proj.get(dir):
                    proj[dir].append(p)
                else:
                    proj[dir] = [p]
    for p1 in dirs:
        for p2 in dirs:
            if p1 != p2:
                Gproj.add_edge(p1,p2)
                        
projects = list(proj.keys())

for p in projects:
    ps = proj[p]
    for p1 in ps:
        for p2 in ps:
            if p1 != p2:
                Gprog.add_edge(p1,p2)

output_rank(Gprog,'prog.txt')
output_rank(Gproj,'proj.txt')
