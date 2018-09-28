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
import os, string, re, glob, time, types, sys

import SCons

# The following adds all SCons SConscript API to the globals of this module.
version = list(map(int,string.split(SCons.__version__,'.')[:3]))
if version[0] >= 1 or version[1] >= 97 or (version[1] == 96 and version[2] >= 90):
    from SCons.Script import *
else:
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

import rsf.path, rsf.sftour, rsf.tex, rsf.conf, rsf.prog

#############################################################################
# CUSTOM BUILDERS
#############################################################################

command = r'\{(?:([^\}\{]*)(?:[^\{\}]*\{[^\}]*\}[^\{\}]*)*)\}'
# matches the contents of {} allowing a single level of nesting
# If an argument of has more than one level of nesting,
# this will not work (Joel)

re_author = re.compile(r'\\author\s*'+command)
re_title = re.compile(r'\\(?:title|chapter)\s*(?:\[[^\]]+\])?\s*'+command)
re_bioname = re.compile(r'{\\bf\s*([^\\\}]+)')

report = string.upper(os.path.basename(os.getcwd()))
group = re.sub('\d+$','',report)

misc = {'pub.tex': 'Articles published or in press',
        'phone.tex': '%s phone directory' % group,
        'bio.tex': 'Research personnel',
        'spons.tex': 'Sponsors'}

full = {
    'SEP': 'Stanford Exploration Project'
    }

prefs = ('abs','pref','ack')

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
        author = re.sub(r'\s+and(?:\s+and)*\s+',' and ',author)
        author = re.sub(r'^\s+','',author)
        author = re.sub(r'\s+$','',author)
    authors = author.split(' and ')
    if len(authors) < 3:
        return author
    else:
        return ', '.join(authors[:-1]) + ', and ' + authors[-1]

def get_title(source,default,tag):
    title = default.get(tag[1])
    if not title:
        paper = source.get_contents()
        # remove comments
        paper = re.sub(r'[%][^\n]+','',paper)
        title = re_title.search(paper)
        if title:
            title = re.sub(r'\\',' ',title.group(1)) # remove line breaks 
    return title

def get_authors(source,default):
    authors = {}
    for src in source:
        tag = paper_tag(str(src))
        author = get_author(src,default,tag)
        if author:
            print("%s: '%s'" % (tag[1],author))
            for person in re.split(r'\s*(?:[\,]|[\,]?\s*\band\b)\s*',author):
                names = string.split(person)
                if names:
                    person = names.pop(0) # first name
                    if names:
                        last = names.pop() # last name
                        person = string.join((person,last),'~')
                        authors[person]=last.capitalize()
    all = [(authors[k],k) for k in list(authors.keys())]
    all.sort()
    return all

def list_authors(all):
    lastone = all.pop()
    if len(all) == 0:
        author = lastone[1]
        print("The author is " + author)
    elif len(all) == 1:
        firstone = all.pop()
        author = '%s and %s' % (firstone[1],lastone[1])
        print("The authors are " + author)
    else:                
        author = string.join([k[1] for k in all],', ')
        author = '%s, and %s' % (author,lastone[1])
        print("The authors are " + author)
    return author
    

def report_toc(target=None,source=None,env=None):
    "Build a table of contents from a collection of papers"
    sections = env.get('sections',{})
    authors = env.get('authors',{})
    titles = env.get('titles',{})
    book = env.get('book')
    if book:
        maketitle = ''
    else:
        maketitle = '\\maketitle'
    toc = open(str(target[0]),'w')
    toc.write(string.join(
        ['%% This file is automatically generated, DO NOT EDIT',
         '\\cleardoublepage',
         '\\renewcommand{\REPORT}{%s}' % report, 
         '\\maintitle{\REPORT\ --- TABLE OF CONTENTS}',
         maketitle,
         '%% start entries\n'],'\n'))
    if book:
        author = '\ '
    for src in source:
        tag = paper_tag(str(src))
        paper = src.get_contents()
        # remove comments
        paper = re.sub(r'[%][^\n]+','',paper)
        if not book:
            author = get_author(src,authors,tag)
        title = get_title(src,titles,tag)
        if tag[1] in sections:
            toc.write('\n\\geosection*{%s}\n' % sections[tag[1]])
        if author and title:
            toc.write('\TOCentry[%s]{%s}{\pageref{%s.start}}\n' %
                      (author,title,tag[1]))
        else:
            print("Could not find author or title")
    toc.write('\n\\geosection*{\ }\n')
    year = get_year(env.get('year'))
    misc['pub.tex'] = '%s article published or in press, %s' % (group,year),
    misc['spons.tex'] = '%s sponsors for %s' % (group,year)
    list(map(lambda x:
        toc.write('\TOCentry{%s}{\pageref{%s.start}}\n' %
                  (misc[x],os.path.splitext(x)[0])),
        list(filter(os.path.isfile,list(misc.keys())))))
    toc.write('\n\\cleardoublepage\n')
    toc.close()
    return 0

def report_tpg(target=None,source=None,env=None):
    "Build the title page"
    tpg = open(str(target[0]),'w')
    tpg.write('%% This file is automatically generated, DO NOT EDIT\n\n')
    tpg.write('\\maintitle{%s}\n' % \
              string.upper(env.get('group',full.get(group))))
    title = env.get('title1')
    if title:
        tpg.write('\\vfill\n\title{%s}\n' % title)
    authors = env.get('authors',{})
    tpg.write('\\mainauthor{%s}\n' % list_authors(get_authors(source,authors)))
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
    tpg.close()
    return 0

def thesis_tpg(target=None,source=None,env=None):
    "Build the title page for a thesis"
    tpg = open(str(target[0]),'w')
    tpg.write('%% This file is automatically generated, DO NOT EDIT\n\n')
    author = env.get('author')
    if author:
        tpg.write('\\vfill\n\\author{%s}\n' % author)
    title = env.get('title')
    if title:
        tpg.write('\\vfill\n\\title{%s}\n' % title)
    univ = env.get('univ')
    if univ and univ == 'UT':
        address = env.get('address')
        if address:
            tpg.write('\\address{%s}\n' % address)
    else:
        year = get_year(env.get('year'))
        tpg.write('\\copyrightyear{%s}\n' % year)
        tpg.write('\\dept{Geophysics}\n')
    tpg.close()
    return 0

def book_tpg(target=None,source=None,env=None):
    "Build the title page for a book"
    tpg = open(str(target[0]),'w')
    tpg.write('%% This file is automatically generated, DO NOT EDIT\n\n')
    author = env.get('author')
    if author:
        tpg.write('\\vfill\n\\author{%s}\n' % author)
    title = env.get('title')
    if title:
        tpg.write('\\vfill\n\\title{%s}\n' % title)
    fig = env.get('fig')
    if fig:
        tpg.write('\\renewcommand{\\plotdir}{%s/Fig}\n'
                  '\\vfill\n\\begin{center}\n'
                  '\\plotbox{%s}{%s}\n\\end{center}\n' % tuple(fig))
    year = get_year(env.get('year'))
    tpg.write('\\date{\\copyright\\ \\today}\n')
    tpg.write('\\maketitle\n')
    tpg.write('\\tableofcontents\n')
    tpg.close()
    return 0

def thesis_intro(target=None,source=None,env=None):
    "Build the introductory part of a thesis"
    intro = open(str(target[0]),'w')
    intro.write('%% This file is automatically generated, DO NOT EDIT\n\n')
    univ = env.get('univ')
    committee = env.get('committee')
    if committee:
        if univ == 'UT':
            super = env.get('supervisor')
            if super:
                if type(super) is list:
                    print(super)
                    intro.write('\\supervisor[%s]{%s}\n' % (super[0],super[1]))
                else:
                    intro.write('\\supervisor{%s}\n' % super)
            intro.write('\\committeemembers')
        for i in range(4):
            if len(committee) > i:
                if univ == 'UT':
                    intro.write(('{%s}\n','[%s]\n')[i+1< len(committee)] % committee[i])
                else:
                    intro.write('\\%s{%s}\n' %
                                (('principaladviser','firstreader','secondreader','thirdreader')[i],
                                 committee[i]))
    if univ == 'UT':
        degrees = env.get('degrees')
        if degrees:
            intro.write('\\previousdegrees{%s}\n' % degrees)
	master = env.get('master')
        if master:
            intro.write('\\degree{Master of Science in Geological Sciences}\n\\degreeabbr{M.S.Geo.Sci.}\n')
            intro.write('\\masterthesis\n')
    else:
        intro.write('\\beforepreface\n')
        intro.write('\\newpage \\ \n')

    if univ == 'UT':
        intro.write('\\copyrightpage\n\\commcertpage\n\\titlepage\n')
        dedication=env.get('dedication')
        if dedication:
            intro.write(dedication.join(['\\begin{dedication}\n',
                                         '\n\\end{dedication}\n']))
        if os.path.isfile('ack.tex'):
            print("Found ack.tex")
            intro.write('\\input{ack}\n'.join(['\\begin{acknowledgments}\n',
                                               '\\end{acknowledgments}\n']))
        if os.path.isfile('abs.tex'):
            print("Found abs.tex")
            intro.write('\\utabstract\n\\indent\n\\input{abs}\n')

        intro.write('\\tableofcontents\n\\listoftables\n\\listoffigures\n')   
    else:
        for pref in prefs:
            tex = pref+'.tex'
            if os.path.isfile(tex): # input if file exists
                print("Found %s" % tex)
                intro.write('\\input{%s}\n' % pref)
        intro.write('\\afterpreface\n')
    intro.close()
    return 0

def report_bio(target=None,source=None,env=None):
    "Build author biographies"
    bio = open(str(target[0]),'w')
    bio.write('\\maintitle{Research Personnel}\n')
    bio.write('\\begin{bios}\n')
    bio.write('\\label{bio.start}\n')
    bios = env.get('bios')
    authors = env.get('authors',{})
    authors = get_authors(source,authors)
    for author in authors:
        name = author[1]
        biofile = bios.get(name)
        if biofile:
            pdf = os.path.splitext(biofile)[0]+'.pdf'
            if os.path.isfile(pdf):
                bio.write('\\item\\bioplot{%s}{%%\n' % pdf)
                thisbio = open(biofile)
                bio.writelines(thisbio.readlines())
                thisbio.close()
                bio.write('}\n')
            else:
                print('No picture file for ' + name)
        else:
            print('No biography file for ' + name)
    bio.write('\\end{bios}\n')
    bio.close()
    return 0

def include(file,sep=''):
    tex = file+'.tex'
    if os.path.isfile(tex):
        print("Found %s" % tex)
        return '\\include{%s}\t%s\t\\newpage\n' % (file,sep)
    else:
        return ''

def report_all(target=None,source=None,env=None):
    "Build the main paper"
    all = open(str(target[0]),'w')
    grp = env.get('group',full.get(group))
    rep = env.get('report',report)
    list(map(all.write,
        ['%% This file is automatically generated, DO NOT EDIT\n\n',
         '\\renewcommand{\\REPORT}{%s}\n' % rep,
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
         '\\GEOheader{\\GROUP, \\REPORT, \\today}\n'
         ]))
    all.write('%% start of paper list\n')
    resdirs = env.get('resdirs',{})
    for src in source:
        paper = str(src)
        tag = paper_tag(paper)
        stem = os.path.splitext(paper)[0]
        resdir = resdirs.get(tag[0],'Fig')
        all.write('\\setfigdir{%s}' % resdir)
        all.write('\\GEOpaper{%s}{%s}\t\\include{%s}\n' % (tag[0],tag[1],stem))
        all.write('\\cleardoublepage\n')
    all.write('%% end of paper list\n')
    for tex in list(misc.keys()):
        all.write(include(os.path.splitext(tex)[0]))
    
    index = env.get('index')
    if index:
        list(map(all.write,
            ['\\cleardoublepage\n',
             '\\addcontentsline{toc}{chapter}{Index}\n',
             '\\index{index}\n',
             '\\printindex\n']))
    biblio = env.get('biblio')
    if biblio:
        all.write('\\cleardoublepage\n\\addcontentsline{toc}{chapter}{Bibliography}\n\\bibliographystyle{seg}\n\\bibliography{%s}\n'
                  % biblio)
    if os.path.isfile('vita.tex'):
        print("Found vita.tex")
        all.write('\\begin{vita}\n\\input{./vita}\n\\end{vita}\n')
    all.close()
    return 0

# ------------------------------------------------------------
def parts(target=None,source=None,env=None):
    "Build a book from parts"

    all = open(str(target[0]),'w')
    grp = env.get('group',full.get(group))
    rep = env.get('report',report)

    list(map(all.write,
        ['%% This file is automatically generated, DO NOT EDIT\n\n',
         '\\renewcommand{\\REPORT}{%s}\n' % rep,
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
         '\\GEOheader{\\GROUP, \\REPORT, \\today}\n'
         ]))
    all.write('%% start of paper list\n')
    resdirs = env.get('resdirs',{})

    parts = env.get('parts', None)
    chpts = env.get('chapters', None)
#    for (part,chapters) in zip(parts,chpts):
#        print 'part=', part, ', chapters=',chapters

    nprt = len(parts)
    for iprt in range(nprt):
    
        prt = str(parts[iprt])
        all.write('''\\part{%s}\n\n''' % prt)
 
        for chp in chpts[iprt]:
            paper = str(chp)
            tag = paper_tag(paper)
            stem = os.path.splitext(paper)[0]
            resdir = resdirs.get(tag[0],'Fig')
    
            all.write('\\setfigdir{%s}\n' % resdir)
            all.write('\\GEOpaper{%s}{%s}\n\\include{%s}\n' % (tag[0],tag[1],stem))
            all.write('\\cleardoublepage\n\n')
    all.write('%% end of paper list\n')

    all.close()
    return 0
# ------------------------------------------------------------

def tour(target=None,source=None,env=None):
    dirs = list(map(os.path.dirname,list(map(str,source))))
    command = env.get('command')
    rsf.sftour.tour(dirs,command,0)
    return 0

Tour = Builder(action = Action(tour),varlist=['command'])

def Sections(report):
    sections = {}
    for section in report:
        first = Split(section[1])[0]
        sections[paper_tag(first)[1]]= section[0]
    return sections

class RSFReport(Environment):
    def __init__(self,**kw):
        Environment.__init__(*(self,), **kw)
        rsf.conf.set_options(self)
        
        self.Append(BUILDERS={'Tour':Tour})
        
        cwd = os.getcwd()
        # create a hierarcical structure
        self.tree = (os.path.basename(os.path.dirname(cwd)),
                     os.path.basename(cwd))

        root = rsf.prog.RSFROOT
        self.doc = os.environ.get('RSFBOOK',os.path.join(root,'share','madagascar'))
        for level in self.tree:
            if level:
                self.doc = os.path.join(self.doc,level)
        rsf.path.mkdir(self.doc)
        self.paper = 1
        self.collection = 0
    def Papers(self,papers,**kw):
        self.collection = 1
        # get list of papers
        if type(papers[0]) is tuple:
            sections = Sections(papers)
            kw.update({'sections':sections})
            papers = Split(string.join([x[1] for x in papers]))
        for i in range(len(papers)):
            paper = papers[i]
            if paper[-4:] != '.tex':
                papers[i] = paper + '/paper.tex'
            elif '/' in paper:
                self.paper = 0
        # make table of contents
        kw.update({'action':Action(report_toc),
                   'varlist':['year','sections','authors','titles']})
        self.Command(*('toc.tex',papers), **kw)
        rsf.tex.Paper('toc',lclass='georeport',scons=0)
        list(map(lambda tex: self.Depends('toc.tex',tex),
            list(filter(os.path.isfile,list(misc.keys())))))
        # make title page
        kw.update({'action':Action(report_tpg),
                   'varlist':['group','title1','authors','title2','line',
                              'fig','year','copyr']})
        self.Command(*('tpg.tex',papers), **kw)
        rsf.tex.Paper('tpg',lclass='georeport',scons=0)
        # make biographies
        bios = kw.get('bios')
        biofiles = {}
        if bios:
            for biofile in glob.glob(os.path.join(bios,'*.bio')):
                bio = open(biofile)
                line = bio.readline()
                name = re_bioname.search(line)
                if name:
                    name = string.split(name.group(1))
                    biofiles[string.join([name[0],name[-1]],'~')] = biofile
                bio.close()
        if biofiles:
            kw.update({'action':Action(report_bio),
                       'varlist':['bios'],
                       'bios':biofiles})
            self.Command(*('bio.tex',papers), **kw)
            self.Depends('bio.tex','tpg.tex')
            rsf.tex.Paper('bio',lclass='georeport',scons=0)
        # make report
        kw.update({'action':Action(report_all),'varlist':['group','resdirs']})
        self.papers = papers
        self.Command(*('book.tex',papers), **kw)
        self.Depends('book.tex','toc.tex')
        self.Depends('book.tex','tpg.tex')
        if biofiles:
            self.Depends('book.tex','bio.tex')
        # Touring individual papers
        for target in ('html','pdf','install'):
            if self.paper:
                comm = 'scons -Q ' + target
            else:
                comm = ['scons -Q %s.%s' % \
                           (os.path.splitext(os.path.basename(x))[0],target) for x in papers]
            self.Tour(target+'s',papers,command=comm)
    def Thesis(self,chapters,**kw):
        # get list of chapters
        for i in range(len(chapters)):
            chapter = chapters[i]
            # add suffix if missing
            if chapter[-4:] != '.tex':
                chapters[i] = chapter + '/paper.tex'
        # make title page
        kw.update({'action':Action(thesis_tpg),
                   'varlist':['author','title','year','univ','address']})
        self.Command(*('tpg.tex',None), **kw)
        rsf.tex.Paper('tpg',lclass='stanford-thesis',scons=0)
        # make introductory materials
        kw.update({'action':Action(thesis_intro),
                   'varlist':['supervisor','committee','master',
                              'univ','degrees','dedication']})
        self.Command(*('intro.tex',None), **kw)
        rsf.tex.Paper('intro',lclass='stanford-thesis',scons=0)
        for pref in prefs:
            tex = pref+'.tex'
            if os.path.isfile(tex):
                self.Depends('intro.tex',tex)
        # make report
        kw.update({'action':Action(report_all),
                   'varlist':['group','resdirs','biblio','report']})
        self.papers = chapters
        self.Command(*('book.tex',chapters), **kw)
        self.Depends('book.tex','tpg.tex')
        self.Depends('book.tex','intro.tex')
        for target in ('html','pdf','install'):	
            comm = 'scons -Q ' + target
            self.Tour(target+'s',chapters,command=comm)
    def Book(self,chapters,**kw):
        self.collection = 1 # for now
        # get list of chapters
        for i in range(len(chapters)):
            chapter = chapters[i]
            # add suffix if missing
            if chapter[-4:] != '.tex':
                chapters[i] = chapter + '/paper.tex'
        # make table of contents
        kw.update({'action':Action(report_toc),
                   'varlist':['year','sections','book'],
                   'book':1})
        self.Command(*('toc.tex',chapters), **kw)
        rsf.tex.Paper('toc',lclass='georeport',scons=0)
        # make title page
        kw.update({'action':Action(book_tpg),
                   'varlist':['author','title','fig']})
        self.Command(*('tpg.tex',None), **kw)
        rsf.tex.Paper('tpg',lclass='georeport',options='book',scons=0)
        # make book
        kw.update({'action':Action(report_all),'varlist':['index']})
        self.papers = chapters
        self.Command(*('book.tex',chapters), **kw)
        self.Depends('book.tex','tpg.tex')
        self.Depends('book.tex','intro.tex')
        for target in ('html','pdf','install'):	
            comm = 'scons -Q ' + target
            self.Tour(target+'s',chapters,command=comm)

    # ------------------------------------------------------------
    def PartBook(self,chapters,**kw):
        self.collection = 1 # for now
        print("BOOK")

        # make title page
        kw.update({'action':Action(book_tpg),'varlist':['author','title']})
        self.Command(*('tpg.tex',None), **kw)
        rsf.tex.Paper('tpg',lclass='georeport',options='book',scons=0)

        # ------------------------------------------------------------
        nprt=len(chapters)
        p=[]
        c=[]
        # loop over parts
        for iprt in range(nprt):
            prt =       chapters[iprt][0]  # separate parts
            chp = Split(chapters[iprt][1]) # separate chapters

            # add suffix if missing
            for ichp in range(len(chp)):
                chapter = chp[ichp]
                if chapter[-4:] != '.tex':
                    chp[ichp] = chapter + '/paper.tex'

            p.append(prt) # collect parts
            c.append(chp) # collect chapters (w/ file names)

        # make book
        kw.update({'action':Action(parts), 'parts':p, 'chapters':c})
        self.Command('book.tex', [], **kw)

        self.papers = c

        self.Depends('book.tex','tpg.tex')
        self.Depends('book.tex','intro.tex')

        for target in ('html','pdf','install'):
            comm = 'scons -Q ' + target
            self.Tour(target+'s',c,command=comm)
    # ------------------------------------------------------------

    def End(self,lclass='georeport',**kw):
        kw.update({'lclass':lclass})
        rsf.tex.Paper(*('book',), **kw)
        self.Alias('pdf','book.pdf')
        self.Depends('book.pdf','pdfs')
        self.Depends('book.pdf',self.papers)
        self.Alias('read','book.read')
        self.Alias('print','book.print')
        if self.collection:
            self.Alias('html','toc_html')
            self.Depends('toc_html/index.html','htmls')
            self.Depends('html','htmls')
            html = os.path.join(self.doc,'index.html')
            self.Command(html,'toc_html/index.html',
                         'cd $SOURCE.dir && cp -R * $TARGET.dir && cd ..')
        else:
            self.Alias('html','book_html')
            html = os.path.join(self.doc,'index.html')
            self.Command(html,'book_html/index.html',
                         'cd $SOURCE.dir && cp -R * $TARGET.dir && cd ..')
        self.Install(self.doc,'book.pdf')
        self.Alias('www',self.doc)
        self.Depends('www','installs')
        self.Default('pdf')
        self.Default('book.tex')

# Default report
book = RSFReport()
def Papers(papers=glob.glob('[a-z]*/paper.tex'),**kw):
    return book.Papers(*(papers,), **kw)
def Thesis(chapters=glob.glob('[a-z]*/paper.tex'),**kw):
    return book.Thesis(*(chapters,), **kw)
def Book(chapters=glob.glob('[a-z]*/paper.tex'),**kw):
    return book.Book(*(chapters,), **kw)
def PartBook(chapters=glob.glob('[a-z]*/paper.tex'),**kw):
    return book.PartBook(*(chapters,), **kw)
def End(lclass='georeport',**kw):
    return book.End(*(lclass,), **kw)



