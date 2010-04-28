#!/usr/bin/env python
'''
Establish comparative ranking of different programs and projects
'''

from rsfdoc import progs
import rsfprog

import networkx as nx

def output_rank(G,filename):
    rank = nx.pagerank(G)
    keys = rank.keys()
    keys.sort(lambda x,y: cmp(rank[y],rank[x]))

    ofile = open(filename,'w')
    for key in keys:
        ofile.write('%s %s\n' % (rank[key],key))
    ofile.close()
    

Gprog = nx.Graph()
Gproj = nx.Graph()

programs = progs.keys()

proj = {}
for p in programs:
    uses = progs[p].uses
    list = []
    for book in uses.keys():
        for chapter in uses[book].keys():
            for project in uses[book][chapter]:
                dir = '/'.join([book,chapter,project])
                list.append(dir)

                if proj.get(dir):
                    proj[dir].append(p)
                else:
                    proj[dir] = [p]

                for project2 in uses[book][chapter]:
                    if project2 != project:
                        dir2 = '/'.join([book,chapter,project2])
                        Gproj.add_edge(dir,dir2)
                        
projects = proj.keys()

for p in projects:
    ps = proj[p]
    for p1 in ps:
        for p2 in ps:
            if p1 != p2:
                Gprog.add_edge(p1,p2)

output_rank(Gprog,'prog.txt')
output_rank(Gproj,'proj.txt')





