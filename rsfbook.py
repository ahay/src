##   Copyright (C) 2004 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
import os, string, re, glob

# The following adds all SCons SConscript API to the globals of this module.
import SCons.Script.SConscript
globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

import rsftex

#############################################################################
# CUSTOM BUILDERS
#############################################################################

command = r'\{(?:([^\}\{]*)(?:[^\{\}]*\{[^\}]*\}[^\{\}]*)*)\}'
# matches the contents of {} allowing a single level of nesting
# If an argument of has more than one level of nesting,
# this will not work (Joel)

re_author = re.compile(r'\\author\s*'+command)
r_title = re.compile(r'\\(?:title|chapter)\s*(?:\[[^\]]+\])?\s*'+command)

def report_toc(target=None,source=None,env=None):
    "Build a table of contents from a collection of papers"
    toc = open(str(target[0]),'w')
    report = string.replace(string.upper(os.path.basename(os.getcwd())),
                            'SEP','SEP--')
    toc.write('''
    \cleardoublepage
    \renewcommand{\REPORT}{%s} 
    \title{\REPORT\ --- TABLE OF CONTENTS}
    \maketitle
    %% start entries

    ''' % report)
    ### Do sections later
    for src in source:
        print "paper=" + str(src)
        paper = src.get_contents()
        author = re_author.search(paper)
        if author:
            print "author=" + author.group(1)
    toc.close()
    return 0

Toc = Builder(action = Action(report_toc))

class Report(Environment):
    def __init__(self,papers=glob.glob('[a-z]*/paper.tex'),**kw):
        apply(Environment.__init__,(self,),kw)
        self.Append(BUILDERS={'Toc':Toc})
        self.Toc('toc.tex',papers)

default = Report()

