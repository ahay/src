#!/usr/bin/env python

# Copyright (C) 2019 University of Texas at Austin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import nbformat, sys, re

scons = {}

def ipynb(figs,notebook):
    'write a notebook'
    nb2 = nbformat.v4.new_notebook()
    nb2.cells.append(nbformat.v4.new_code_cell('from m8r import view'))
    n = 0
    for fig in figs:
        n += 1
        code = '%%%%file %d_%s.scons \n\n%s' % (n,fig,scons[fig])           
        nb2.cells.append(nbformat.v4.new_code_cell(code))
        if fig != 'whatever':
            nb2.cells.append(nbformat.v4.new_code_cell('view("%s")'  % fig))
    nbformat.write(nb2,notebook)

def parse(sconstruct):
    'parse SConstruct, extract code for each figure'
    result = re.compile('Result\s*\(\s*[\'\"](?P<fig>[^\'\"]+)')
    figs = []
    block = ''
    brackets = 0
    for line in sconstruct:
        if brackets: # unclosed Result()
            for c in line:
                if c == '(':
                    brackets += 1
                elif c == ')':
                    brackets -= 1
            if brackets == 0:
                block += line.rstrip()
                scons[fig] = block
                block = ''
            else:
                block += line
            continue
        
        if re.search('from rsf.proj import \*\s*',line):
            continue
        if re.search('End()\s*',line):
            continue
        
        res = result.search(line)
        if res:
            fig = res.group('fig')
            figs.append(fig)
            
            for c in line:
                if c == '(':
                    brackets += 1
                elif c == ')':
                    brackets -= 1
            if brackets == 0:
                block += line.rstrip()
                scons[fig] = block
                block = ''
            else:
                block += line
        else:
            block += line
    block = block.rstrip()
    if block:
        fig = 'whatever'
        figs.append(fig)
        scons[fig] = block
    return figs
    
if __name__ == "__main__":
    figs = parse(sys.stdin)
    ipynb(figs,sys.stdout)
    sys.exit(0)
