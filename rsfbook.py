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
import os, string, re, glob, time

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
re_title = re.compile(r'\\(?:title|chapter)\s*(?:\[[^\]]+\])?\s*'+command)

report = string.upper(os.path.basename(os.getcwd()))
group = re.sub('\d+$','',report)

misc = {'pub.tex': 'Articles published or in press',
        'phone.tex': '%s phone directory' % group,
        'bio.tex': 'Research personnel',
        'spons.tex': 'Sponsors'}

def report_toc(target=None,source=None,env=None):
    "Build a table of contents from a collection of papers"
    toc = open(str(target[0]),'w')
    toc.write(string.join(
        ['%% This file is automatocally generated, DO NOT EDIT',
         '\\cleardoublepage',
         '\\renewcommand{\REPORT}{%s}' % report, 
         '\\title{\REPORT\ --- TABLE OF CONTENTS}',
         '\\maketitle',
         '%% start entries\n'],'\n'))
    ### Do sections later
    for src in source:
        dir = os.path.basename(os.path.dirname(str(src)))
        paper = src.get_contents()
        author = re_author.search(paper)
        title = re_title.search(paper)
        if author and title:
            title = re.sub(r'\\',' ',title.group(1)) # remove line breaks 
            toc.write('\TOCentry[%s]{%s}{\pageref{%s.start}}\n' %
                      (author.group(1),title,dir))
        else:
            print "Could not find author or title"
    year = env.get('year')
    if not year:
        year = time.localtime(time.time())[0]
    year = str(year)
    year = string.join([str(int(year)-1),year[-2:]],'-')
    misc['pub.tex'] = '%s article published or in press, %s' % (group,year),
    misc['spons.tex'] = '%s sponsors for %s' % (group,year)
    map(lambda x:
        toc.write('\TOCentry{%s}{\pageref{%s.start}}\n' %
                  (misc[x],os.path.splitext(x)[0])),
        filter(os.path.isfile,misc.keys()))
    toc.write('\n\\cleardoublepage\n')
    toc.close()
    return 0

Toc = Builder(action = Action(report_toc), varlist=['year'])

class RSFReport(Environment):
    def __init__(self,**kw):
        apply(Environment.__init__,(self,),kw)
        self.Append(BUILDERS={'Toc':Toc})        
    def Papers(self,papers=glob.glob('[a-z]*/paper.tex'),year=None):
        self.Toc('toc.tex',papers,year=year)
        map(lambda tex: self.Depends('toc.tex',tex),
            filter(os.path.isfile,misc.keys()))

# Default report
book = RSFReport()
def Papers(papers=glob.glob('[a-z]*/paper.tex'),**kw):
    return apply(book.Papers,(papers,),kw)
    
