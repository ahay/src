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
import os, string, re, glob, time, types

import SCons

# The following adds all SCons SConscript API to the globals of this module.
version = map(int,string.split(SCons.__version__,'.'))
if version[1] >= 97 or (version[1] == 96 and version[2] >= 90):
    from SCons.Script import *
else:
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

import rsftex, sftour

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

full = {
    'SEP': 'Stanford Exploration Project'
    }

def get_year(default):
    year = default
    if not year:
        year = time.localtime(time.time())[0]
        year = str(year)
        year = string.join([str(int(year)-1),year[-2:]],'-')
    return year

def paper_tag(paper):
    return (os.path.basename(os.path.dirname(paper)),
            re.sub('-paper$','',
                   re.sub(r'/([^/]+)\.tex$',r'-\1',paper)))

def get_author(source,default,tag):
    author = default.get(tag[1])
    if not author:
        paper = source.get_contents()
        # remove comments
        paper = re.sub(r'[%][^\n]+','',paper)
        author = re_author.search(paper)
        if author:
            author = author.group(1)
    if author:
        author = re.sub(r'(?:\,|\;|\\\&)',' and ',author)
        author = re.sub(r'\s+and\s+and\s+',' and ',author)
        author = re.sub(r'^\s+','',author)
        author = re.sub(r'\s+$','',author)
    return author

def get_authors(source,default):
    authors = {}
    for src in source:
        tag = paper_tag(str(src))
        author = get_author(src,default,tag)
        if author:
            print "%s: '%s'" % (tag[1],author)
            for person in re.split(r'\s*\band\b\s*',author):
                names = string.split(person)
                if names:
                    person = names.pop(0) # first name
                    if names:
                        last = names.pop() # last name
                        person = string.join((person,last),'~')
                        authors[person]=last
    all = map(lambda k: (authors[k],k),authors.keys())
    all.sort()
    lastone = all.pop()
    if len(all) == 0:
        author = lastone[1]
        print "The author is " + author
    elif len(all) == 1:
        firstone = all.pop()
        author = '%s and %s' % (firstone[1],lastone[1])
        print "The authors are " + author
    else:                
        author = string.join(map(lambda k: k[1],all),', ')
        author = '%s, and %s' % (author,lastone[1])
        print "The authors are " + author
    return author

def report_toc(target=None,source=None,env=None):
    "Build a table of contents from a collection of papers"
    toc = open(str(target[0]),'w')
    toc.write(string.join(
        ['%% This file is automatocally generated, DO NOT EDIT',
         '\\cleardoublepage',
         '\\renewcommand{\REPORT}{%s}' % report, 
         '\\maintitle{\REPORT\ --- TABLE OF CONTENTS}',
         '\\maketitle',
         '%% start entries\n'],'\n'))
    sections = env.get('sections',{})
    authors = env.get('authors',{})
    for src in source:
        tag = paper_tag(str(src))
        paper = src.get_contents()
        # remove comments
        paper = re.sub(r'[%][^\n]+','',paper)
        author = get_author(src,authors,tag)
        title = re_title.search(paper)
        if sections.has_key(tag[1]):
            toc.write('\n\\geosection*{%s}\n' % sections[tag[1]])
        if author and title:
            title = re.sub(r'\\',' ',title.group(1)) # remove line breaks 
            toc.write('\TOCentry[%s]{%s}{\pageref{%s.start}}\n' %
                      (author,title,tag[1]))
        else:
            print "Could not find author or title"
    year = get_year(env.get('year'))
    misc['pub.tex'] = '%s article published or in press, %s' % (group,year),
    misc['spons.tex'] = '%s sponsors for %s' % (group,year)
    map(lambda x:
        toc.write('\TOCentry{%s}{\pageref{%s.start}}\n' %
                  (misc[x],os.path.splitext(x)[0])),
        filter(os.path.isfile,misc.keys()))
    toc.write('\n\\cleardoublepage\n')
    toc.close()
    return 0

def report_tpg(target=None,source=None,env=None):
    "Build the title page"
    tpg = open(str(target[0]),'w')
    tpg.write('%% This file is automatocally generated, DO NOT EDIT\n\n')
    tpg.write('\\maintitle{%s}\n' % \
              string.upper(env.get('group',full.get(group))))
    title = env.get('title1')
    if title:
        tpg.write('\\vfill\n\title{%s}\n' % title)
    authors = env.get('authors',{})
    tpg.write('\\mainauthor{%s}\n' % get_authors(source,authors))
    title = env.get('title2')
    if title:
        tpg.write('\\vfill\n\maintitle{%s}\n' % title)
    line = env.get('line')
    if line:
        tpg.write('\\vfill\n\\begin{center}\n'
                  '\\bfseries%%\n%s\n\\end{center}\n' % line)
    fig = env.get('fig')
    if fig:
        tpg.write('\\renewcommand{\\plotdir}{%s/Fig}\n'
                  '\\vfill\n\\begin{center}\n'
                  '\\plotbox{%s}{%s}\n\\end{center}\n' % tuple(fig))
    year = get_year(env.get('year'))
    copyr = env.get('copyr')
    tpg.write('\n\\newpage\\GEOcopyr{%s}{%s}\n' % (year,copyr))
    return 0

def include(file,sep=''):
    tex = file+'.tex'
    if os.path.isfile(tex):
        print "Found %s" % tex
        return '\\include{%s}\t%s\t\\newpage\n' % (file,sep)
    else:
        return ''

def report_all(target=None,source=None,env=None):
    "Build the main paper"
    all = open(str(target[0]),'w')
    grp = env.get('group',full.get(group))
    map(all.write,
        ['%% This file is automatocally generated, DO NOT EDIT\n\n',
         '\\renewcommand{\\REPORT}{%s}\n' % report,
         '\\renewcommand{\\GROUP}{%s}\n' % grp,
         '\\renewcommand{\\thepage}{}\n',
         include('tpg','\\cleardoublepage'),
         include('preface'),
         '\\pagenumbering{roman}\n',
         '\\setcounter{page}{1}\n',
         include('toc'),
         include('intro'),
         '\\cleardoublepage\n',
         '\\pagenumbering{arabic}\n',
         '\\setcounter{page}{1}\n',
         '\\GEOheader{\\GROUP, Report \\REPORT, \\today}\n'
         ])
    all.write('%% start of paper list\n')
    resdirs = env.get('resdirs',{})
    for src in source:
        paper = str(src)
        tag = paper_tag(paper)
        stem = os.path.splitext(paper)[0]
        resdir = resdirs.get(tag[0],'Fig')
        all.write('\\setfigdir{%s}' % resdir)
        all.write('\\GEOpaper{%s}{%s}\t\\include{%s}\n' % (tag[0],tag[1],stem))
        all.write('\\cleardoublepage')
    all.write('%% end of paper list\n')
    map(all.write, [ include('post'), ])
    return 0

def tour(target=None,source=None,env=None):
    dirs = map(os.path.dirname,map(str,source))
    command = env.get('command')
    sftour.tour(dirs,command,0)
    return 0

Tour = Builder(action = Action(tour),varlist=['command'])

def Sections(report):
    sections = {}
    for section in report:
        first = Split(section[1])[0]
        sections[first]= section[0]
    return sections

class RSFReport(Environment):
    def __init__(self,**kw):
        apply(Environment.__init__,(self,),kw)
        self.Append(BUILDERS={'Tour':Tour})
        cwd = os.getcwd()
        # create a hierarcical structure
        self.tree = (os.path.basename(os.path.dirname(cwd)),
                     os.path.basename(cwd))
        self.doc = os.environ.get(
            'RSFDOC',os.path.join(os.environ.get('RSFROOT'),'doc'))
        for level in self.tree:
            if level:
                self.doc = os.path.join(self.doc,level)
        rsftex.mkdir(self.doc)
        self.paper = 1
    def Papers(self,papers,**kw):
        if type(papers[0]) is types.TupleType:
            kw.update({'sections':Sections(papers)})
            papers = Split(string.join(map(lambda x: x[1],papers)))
        for i in range(len(papers)):
            paper = papers[i]
            if paper[-4:] != '.tex':
                papers[i] = paper + '/paper.tex'
            elif '/' in paper:
                self.paper = 0
        kw.update({'action':Action(report_toc),
                   'varlist':['year','sections','authors']})
        apply(self.Command,('toc.tex',papers),kw)
        rsftex.Paper('toc',lclass='georeport',scons=0)
        kw.update({'action':Action(report_tpg),
                   'varlist':['group','title1','authors','title2','line',
                              'fig','year','copyr']})
        apply(self.Command,('tpg.tex',papers),kw)
        rsftex.Paper('tpg',lclass='georeport',scons=0)
        kw.update({'action':Action(report_all),'varlist':['group','resdirs']})
        apply(self.Command,('book.tex',papers),kw)
        for file in ['tpg.tex','toc.tex']:
            map(lambda tex: self.Depends(file,tex),
                filter(os.path.isfile,misc.keys()))
            self.Depends('book.tex',file)
        # Touring individual papers
        for target in ('html','pdf','install'):
            if self.paper:
                comm = 'scons -Q ' + target
            else:
                comm = map(lambda x: 'scons -Q %s.%s' % \
                           (os.path.splitext(os.path.basename(x))[0],target),
                           papers)
            self.Tour(target+'s',papers,command=comm)
    def End(self,**kw):
        kw.update({'lclass':'georeport'})
        apply(rsftex.Paper,('book',),kw)
        self.Alias('pdf','book.pdf')
        self.Depends('book.pdf','pdfs')
        self.Alias('read','book.read')
        self.Alias('print','book.print')
        self.Alias('html','toc_html')
        self.Depends('toc_html/index.html','htmls')
        self.Depends('html','htmls')
        self.Install(self.doc,'book.pdf')
        html = os.path.join(self.doc,'index.html')
        self.Command(html,'toc_html/index.html',
                    'cd $SOURCE.dir && cp -R * $TARGET.dir && cd ..')
        self.Alias('www',self.doc)
        self.Depends('www','installs')
        self.Default('pdf')

# Default report
book = RSFReport()
def Papers(papers=glob.glob('[a-z]*/paper.tex'),**kw):
    return apply(book.Papers,(papers,),kw)
def End(**kw):
    return apply(book.End,(),kw)



