#!/usr/bin/env python
'''
Find modules used in books
'''

import sys
from rsf.doc import progs
import rsf.prog

if len(sys.argv) != 2:
    print(f'Usage: {sys.argv[0]} book1,book2,...')
    sys.exit(0)

# get all programs
programs = list(progs.keys())
programs.sort()

# get book list
books = sys.argv[1].split(',')

def in_books(p, books):
    "if a program is used in one of the books"
    uses = progs[p].uses
    for book in list(uses.keys()):
        if book in books:
            return True
    return False

for p in programs:
    if in_books(p, books):
        print(p)



