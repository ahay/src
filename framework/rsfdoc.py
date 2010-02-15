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

import pydoc
import re, sys, os, glob, signal
import rsfpath

try:
    import datetime
    have_datetime_module = True
except: # Python < 2.3
    have_datetime_module = False

progs = {}
data = {}

# documented programs
docprogs = '''
add attr cat cmplx conjgrad cp cut dd disfil dottest get headerattr
headercut headermath headersort headerwindow in interleave mask math
pad prep4plot put real remap1 reverse rm rotate rtoc scale segyread
segywrite spike spray stack stretch transp window sizes figlist booklist
'''.split()

def handler(signum, frame):
    'signal handler for abortion [Ctrl-C]'
    sys.stderr.write('\n[Ctrl-C] Aborting...\n')
    if child:
        os.kill (signal.SIGINT,child)
    sys.exit(-1)

signal.signal(signal.SIGINT,handler) # handle interrupt

child = None
def syswait(comm):
    'Interruptable system command'
    global child
#    print comm
    child = os.fork()
    if child:
        (pid,exit) = os.waitpid(child,0)
        child = 0
        return exit
    else:
        os.system(comm)
#        print "Done"
        os._exit(0)

def subdirs():
    return filter(lambda x: x[-5:] != '_html',
                  filter(os.path.isdir,glob.glob('[a-z]*')))

def getprogs(target=None,source=None,env=None):
    out = open(str(target[0]),'w')
    dirs = env.get('dirs')
    out.write('import sys, os\n\n')
    out.write('import rsfdoc\n\n')
    for mod in dirs:
        out.write('import sf%s\n' % mod)
    out.write('\nimport vpplot\n\n')
    out.write('''
try:
    import rsfuse
except:
    pass

def selfdoc():
    'Display man page'
    prognm = os.path.basename(sys.argv[0])
    if prognm[0] == 'M' and prognm[-3:] == '.py':
        # User testing code in his local directory
        prognm = 'sf' + prognm[1:-3]
        msg = 'Self-doc may be out of synch, do "scons install" in RSFSRC'
        sys.stderr.write(msg)

    prog = rsfdoc.progs.get(prognm)
    if prog != None:
        prog.document()
    else:
        sys.stderr.write('No installed man page for ' + prognm)
''')
    out.close()

def use(target=None,source=None,env=None):
    out = open(str(target[0]),'w')
    doc = 'import rsfdoc\n'
    cwd = os.getcwd()
    bookdir = env.get('book')
    scons = env.get('scons')
    if os.path.isdir(bookdir):
        os.chdir(bookdir)
        cwtop = os.getcwd()
        for book in subdirs():
            os.chdir(book)
            cwbook = os.getcwd()
            print "%s..." % book
            for chapter in subdirs():
                os.chdir(chapter)
                print "...%s" % chapter

                syswait(scons + ' -s uses data')

                datapath = rsfpath.datapath()
                path = os.path.dirname(datapath)
                if datapath[:2] != './':
                    path = os.path.join(path,book,chapter)
                uses = os.path.join(path,'.sf_uses')

                if os.path.isfile(uses):
                    sout = open(uses,'r')
                    doc = doc + sout.read()
                    sout.close()

                os.chdir(cwbook)
            os.chdir(cwtop)
        os.chdir(cwd)
    out.write(doc)
    out.close()
    return 0

def selfdoc(target=None,source=None,env=None):
    src = str(source[0])
    doc = open(str(target[0]),"w")
    rsfprefix = env.get('rsfprefix','sf')
    rsfsuffix = env.get('rsfsuffix','rsf')
    lang = env.get('lang','c')
    getprog(src,doc,lang,rsfprefix,rsfsuffix)
    doc.close()
    return 0

def bold(text):
    """Format a string in bold by overstriking."""
    return ''.join(map(lambda ch: ch + "\b" + ch, text))

def underline(text):
    """Format a string in underline by overstriking."""
    return ''.join(map(lambda ch: ch + "\b_", text))

def replace(text, *pairs):
    """Do a series of global replacements on a string."""
    while pairs:
        text = pairs[1].join(text.split(pairs[0]))
        pairs = pairs[2:]
    return text

def section(head,body):
    text = "\n".join(map(lambda line: "\t" + line,
                         body.split("\n")))
    return bold(head.upper()) + "\n" + text + "\n"

def page(title, contents):
    """Format an HTML page."""
    return '''
    <!doctype html PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
    <html><head><title>RSF: %s</title>
    <style type="text/css"><!--
    TT { font-family: lucidatypewriter, lucida console, courier }
    --></style></head><body bgcolor="#f0f0f8">
    %s
    </body></html>''' % (title, contents)

def heading(title, fgcol, bgcol,extras='',add='&nbsp;<br>'):
    """Format a page heading."""
    return '''
    <table width="100%%" cellspacing=0 cellpadding=2 border=0
    summary="heading">
    <tr bgcolor="%s">
    <td valign=bottom>%s
    <font color="%s" face="helvetica, arial">&nbsp;<br>%s</font></td
    ><td align=right valign=bottom
    ><font color="%s" face="helvetica, arial">%s</font></td></tr></table>
    ''' % (bgcol, add,fgcol, title, fgcol, extras or '&nbsp;')

def html_section(title, fgcol, bgcol, contents, width=6,
                 prelude='', marginalia=None, gap='&nbsp;'):
    """Format a section with a heading."""
    if marginalia is None:
        marginalia = '<tt>' + '&nbsp;' * width + '</tt>'
    result = '''
    <p>
    <table width="100%%" cellspacing=0 cellpadding=2 border=0 summary="section">
    <tr bgcolor="%s">
    <td colspan=3 valign=bottom>&nbsp;<br>
    <font color="%s" face="helvetica, arial">%s</font></td></tr>
    ''' % (bgcol, fgcol, title)
    if prelude:
        result = result + '''
        <tr bgcolor="%s"><td rowspan=2>%s</td>
        <td colspan=2>%s</td></tr>
        <tr><td>%s</td>''' % (bgcol, marginalia, prelude, gap)
    else:
        result = result + '''
        <tr><td bgcolor="%s">%s</td><td>%s</td>''' % (bgcol, marginalia, gap)
    return result + '\n<td width="100%%">%s</td></tr></table>' % contents

def bigsection(title, *args):
    """Format a section with a big heading."""
    title = '<big><strong>%s</strong></big>' % title
    return apply(html_section,(title,)+args)

def multicolumn(list, format, cols=4):
    """Format a list of items into a multi-column list."""
    result = ''
    rows = (len(list)+cols-1)/cols
    for col in range(cols):
        result = result + '<td width="%d%%" valign=top>' % (100/cols)
        for i in range(rows*col, rows*col+rows):
            if i < len(list):
                result = result + format(list[i]) + '<br>\n'
        result = result + '</td>'
    return '<table width="100%%" summary="list"><tr>%s</tr></table>' % result

class rsfpar(object):
    def __init__(self,type,default='',range='',desc=''):
        self.type = type + ' ' * (len('strings')-len(type))
        if (type == 'bool'):
            if (default == 'true'):
                default = 'y'
            elif (default == 'false'):
                default = 'n'
            range = '[y/n]'
        self.default = "=" + str(default)
        self.range = " " + range
        self.desc = "\t" + desc + "\n"
    def check(self,name,val):
        '''check if val is an acceptable value, return 1 if not'''
        if self.type == 'bool':
            if val != '0' or val != '1' or val != 'y' or val != 'n':
                sys.stderr.write("%s is bool, "
                                 "%s is not acceptable "
                                 "(use [0,1,y,n])\n" % (name,val))
                return 1
        return 0
    def show(self,name):
        return underline(self.type) + " " + \
               bold(name + self.default) + self.range + self.desc
    def html(self,name):
        return self.type + " " + \
               "<strong>" + name + self.default + \
               "</strong>" + self.range
    def text(self,name):
        return ' | '.join([self.type,name,self.default,self.range,self.desc])
    def spec(self,name):
        type = self.type.rstrip()
        range = self.range
        default = self.default
        desc = self.desc.replace('"',"'")

        DescriptionIsLong = desc.count('\n') > 1
        desc = desc.rstrip('\n') # save trailing line feeds for later

        isPort = False
        if type.endswith('s'):
            typesuff = '-list'
            type = type.rstrip('s')
        else:
            typesuff = ''

        if type=='file':
            isPort = True
        elif type=='bool':
            type = 'enum-'+type
            range = ' - '
        elif type=='int':
            range = ' - ' # TO DO: int range detection; 'enum-' type for non-trivial int ranges
        """TO DO: float range detection"""

        """Parse invalid defaults to long description"""
        LHS, default = default.split('=',1)
        try:
            try:
                if type=='enum-bool':
                    if default not in ['y','n']:
                        raise TypeError
                elif '(' in default:
                    raise NameError
                else:
                    DefaultIsValid = isinstance(eval(default),eval(type))
            except (SyntaxError,NameError,TypeError): # exception when default contains expressions
                if default:
                    desc = '%s\n(defaults to: %s)' % (desc,default)
                    DescriptionIsLong = True
                default = ' - '
        finally:
            if len(default)==0:
                default = ' - '

        if DescriptionIsLong:
            desc = desc.replace('\n','\nLDesc:  ') # translate additional info into long description
        if (isPort):
            try:
                data, ext = default.split('.') # assume 'default' carries the filename stated in synopsis
            except:
                data = 'unspecified'
                ext = 'rsf'
            if (desc.find('input')>=0):
                rw = 'r'
            elif (desc.find('output')>=0):
                rw = 'w'
            unspec = {'r':'in', 'w':'out'}
            if (data == unspec[rw]):
                data = '(no hint on content)'
            else:
                data = '(containing %s data)' % data
            line = 'Port:   %s %s %s %s \t%s %s\n' % (name,ext,rw,default,desc,data)
        else:
            type = type + typesuff
            if range.isspace():
                range=' - '
            line = 'Param:  %s %s %s %s \t%s\n' % (name,type,range,default,desc)
        return line
    def mwiki(self,name):
        desc = self.desc.replace('\n','<br>')
        desc = re.sub(r'<br>\s+(\w)',r'\n:\1',desc)
        desc = re.sub(r'(?:\s*<br>)+$','',desc)
        return "|-\n| ''%s'' || '''%s%s''' || %s || %s\n" % \
               (self.type,name,self.default,self.range,desc)
    def latex(self,name):
        tex = '\\underline{%s} & \\textbf{%s%s} & %s & ' % \
              (self.type,name,self.default,self.range)
        return tex + self.desc + '\\\\ \n'
    def man(self,name):
        troff = '.TP\n.I %s\n.B %s\n.B %s\n.R %s%s\n' % \
            (self.type,name,self.default,self.range,self.desc)
        return troff
class rsfdata(object):
    def __init__(self,name):
        self.name = name
        self.uses = {}
    def use (self,book,chapter,project):
        if not self.uses.has_key(book):
            self.uses[book]={}
        if not self.uses[book].has_key(chapter):
            self.uses[book][chapter] = []
        self.uses[book][chapter].append(project)

class rsfprog(object):
    def __init__(self,name,file,desc=None):
        self.name = name
        self.file = file
        self.desc = desc
        self.snps = None
        self.cmts = None
        self.also = None
        self.vers = None
        self.wiki = None
        self.uses = {}
        self.pars = {}
    def check(self,par):
        '''par is name=value pair, check if it makes sense, return 1 if not'''
        pair = par.split('=')
        if len(pair) != 2:
            return 0
        key = pair[0]
        val = pair[1]
        mypar = self.pars.get(key)
        if mypar:
            return mypar.check(key,val)
        if key[-1] in '123456789':
            mypar = self.pars.get(key[:-1]+'#')
            if mypar:
                return mypar.check(key,val)
        sys.stderr.write('No parameter "%s" in %s\n' % (key,self.name))
        return 1
    def weblink(self,wiki):
        self.wiki = wiki
    def version(self,vers):
        self.vers = vers
    def synopsis (self,snps,cmts):
        self.snps = snps
        self.cmts = cmts
    def par (self,name,par):
        self.pars[name] = par
    def use (self,book,chapter,project):
        if not self.uses.has_key(book):
            self.uses[book]={}
        if not self.uses[book].has_key(chapter):
            self.uses[book][chapter] = []
        self.uses[book][chapter].append(project)
    def docstring(self):
        doc = section('name',self.name)
        if self.desc:
            doc = doc + section('description',self.desc)
        if self.snps:
            doc = doc + section('synopsis',self.snps)
        if self.cmts:
            doc = doc + section('comments',self.cmts)
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            pardoc = ''
            for par in pars:
                pardoc = pardoc + self.pars[par].show(par)
            doc = doc + section('parameters',pardoc.rstrip())
        if self.also:
            doc = doc + section('see also',self.also)
        books = self.uses.keys()
        if books:
            usedoc = ''
            books.sort()
            for book in books:
                chapters = self.uses[book].keys()
                chapters.sort()
                for chapter in chapters:
                    for project in self.uses[book][chapter]:
                        usedoc = usedoc + '%s/%s/%s\n' % (book,chapter,project)
            doc = doc + section('used in',usedoc.rstrip())
        doc = doc + section('source',self.file)
        if self.wiki:
            doc = doc + section('documentation',underline(self.wiki))
        if self.vers:
            doc = doc + section('version',self.vers)
        return doc
    def document(self):
        doc = self.docstring()
        pydoc.pager(doc)
    def mwiki(self,dir,name=None):
        if not name:
            name = self.name
        file = open (os.path.join(dir,name + '.wiki'),'w')
        contents = '==%s==\n{| class="wikitable" ' % name
        contents += 'align="center" cellspacing="0" border="1"\n'
        desc = '! colspan="4" style="background:#ffdead;" | %s\n' % self.desc
        contents = contents + desc
        if self.snps:
            contents = contents + '|-\n! colspan="4" | %s\n' % self.snps
        if self.cmts:
            cmts =self.cmts.replace('\n','<br>')
            cmts = re.sub(r'(?:\s*<br>)+$','',cmts)
            contents = contents + '|-\n|  colspan="4" | %s\n' % cmts
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            for par in pars:
                contents = contents + self.pars[par].mwiki(par)
        contents = contents+'|}\n'
        file.write(contents)
        file.close()
    def man(self,dir,name=None):
        if not name:
            name = self.name
        file = open (os.path.join(dir,name + '.1'),'w')
        contents = '.TH %s 1  ' % name
        if have_datetime_module:
            day = datetime.datetime.now()
            month = day.strftime('%B %Y').upper()
            contents += '"%s" ' % month
        contents += 'Madagascar "Madagascar Manuals"\n'
        desc = '.SH NAME\n%s \- %s\n' % (name,self.desc)
        contents = contents + desc
        if self.snps:
            contents = contents + '.SH SYNOPSIS\n.B %s\n' % self.snps
        if self.cmts:
            contents = contents + '.SH COMMENTS\n%s\n' % self.cmts
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            for par in pars:
                contents = contents + self.pars[par].man(par)
        books = self.uses.keys()
        if books:
            usedoc = ''
            books.sort()
            for book in books:
                chapters = self.uses[book].keys()
                chapters.sort()
                for chapter in chapters:
                    for project in self.uses[book][chapter]:
                        usedoc = usedoc + \
                            '.TP\n.I %s/%s/%s\n' % (book,chapter,project)
            if usedoc:
                contents = contents + '.SH USED IN\n%s' % usedoc
        contents = contents + '.SH SOURCE\n.I %s\n' % self.file
        if self.wiki:
            contents = contents + '.SH DOCUMENTATION\n.BR %s\n' % self.wiki
        if self.vers:
            contents = contents + '.SH VERSION\n%s\n' % self.vers
        file.write(contents)
        file.close()
    def latex(self,dir,name=None):
        if not name:
            name = self.name
        file = open (os.path.join(dir,name + '.tex'),'w')
        contents = '\\footnotesize\n'
        name = '\\subsection{%s: %s}\n' % (name,self.desc)
        contents = contents + name
        if self.snps:
            contents = contents + '\\texttt{%s}\n' % self.snps
        if self.cmts:
            contents = contents + '\\begin{verbatim}%s\\end{verbatim}\n' % \
                       self.cmts
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            contents = contents + '\\par\\noindent\n' + \
                       '\\begin{tabular}{p{1in}p{1in}p{1in}p{2.5in}}\n'
            for par in pars:
                contents = contents + self.pars[par].latex(par)
            contents = contents + '\\end{tabular}\n'
        contents = re.sub('([_#])',r'\\\1',contents)
        file.write(contents)
        file.close()
    def text(self,dir,name=None):
        if not name:
            name = self.name
        file = open (os.path.join(dir,name + '.txt'),'w')
        contents = 'Program %s | %s\n' % (name,self.desc)
        if self.snps:
            contents = contents + '[SYNOPSIS]\n%s\n' % self.snps
        if self.cmts:
            contents = contents + '[COMMENTS]\n%s\n' % self.cmts
        pars =  self.pars.keys()
        if pars:
            contents = contents + '[PARAMETERS]\n'
            pars.sort()
            for par in pars:
                contents = contents + self.pars[par].text(par)
        filedir = os.path.split(self.file)[0]
        if filedir:
            contents = contents + '[DIRECTORY]\n%s\n' % filedir
        file.write(contents)
        file.close()
    def spec(self,dir,name=None):
        if not name:
            name = self.name
        file = open (os.path.join(dir,name + '.spec'),'w')
        filedir = os.path.split(self.file)[0]
        doccmd = '%s | cat' % name
        contents = '''[%s]
Cat:    RSF/%s
Desc:   %s
DocCmd: %s
''' % (name,filedir,self.desc.rstrip(' .'),doccmd)
        """ process stdin and stdout ports hidden in synopsis """
        if self.snps:
            tokens = self.snps.split()
            for cue in ['<','>']:
                try:
                    filename = tokens[tokens.index(cue)+1]
                except:
                    continue
                data, ext  = filename.split('.')
                if data in ['in','out']:
                    data = '(no hint on content)'
                else:
                    data = '(containing %s data)' % data
                if cue=='<':
                    line = 'Port:   stdin  %s r req \t%s standard input %s\n' % (ext,ext.upper(),data)
                else:
                    line = 'Port:   stdout %s w req \t%s standard output %s\n' % (ext,ext.upper(),data)
                contents = contents + line
        """TO DO: process comments."""
        pars =  self.pars.keys()
        ParamLines = ''
        if pars:
            pars.sort()
            for par in pars:
                line = self.pars[par].spec(par)
                head, tail = line.split(':',1)
                if head=='Port':
                    contents = contents + line
                else:
                    ParamLines = ParamLines + line
        contents = contents + ParamLines + '\n'
        file.write(contents)
        file.close()
    def html(self,dir,rep):
        file = open (os.path.join(dir,self.name + '.html'),'w')
        name = '<big><big><strong>%s</strong></big></big>' % self.name
        if self.vers:
            name = name + " (%s)" % self.vers
        if self.wiki:
            wiki = '<br><a href="%s">Documentation</a>' % self.wiki
        else:
            wiki = ''
        contents = heading(name,'#ffffff','#7799ee',
                           '<a href="./index.html">index</a><br>'+
                           '<a href="%s/%s?view=markup">%s</a>%s' %
                           (rep,self.file,self.file,wiki))
        if self.desc:
            contents = contents + self.desc
        if self.snps:
            contents = contents + \
                       bigsection('Synopsis','#fffff', '#aa55cc',self.snps)
        if self.cmts:
            contents = contents + self.cmts.replace('\n','<br>\n')
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            pardoc = ''
            bgcol = '#ffc8d8'
            for par in pars:
                pardoc = pardoc + \
                         heading(self.pars[par].html(par),'#000000', bgcol,
                                 self.pars[par].desc.replace('\n','<br>\n'),
                                 add='') + '\n'
                if bgcol=='#ffc8d8':
                    bgcol ='#f0f0f8'
                else:
                    bgcol = '#ffc8d8'
            contents = contents + \
                       bigsection('Parameters','#ffffff', '#ee77aa',pardoc)
        books = self.uses.keys()
        if books:
            usedoc = ''
            books.sort()
            for book in books:
                bookdoc = ''
                chapters = self.uses[book].keys()
                chapters.sort()
                for chapter in chapters:
                    for project in self.uses[book][chapter]:
                        proj = os.path.join(chapter,project)
                        bookdoc = bookdoc + '''
                        <a href="book/%s/%s.html">%s</a><br>
                        ''' % (book,proj,proj)
                usedoc = usedoc + \
                         bigsection(book.upper(),'#000000','#ffd8c8',bookdoc)
            contents = contents + \
                       bigsection('Used In','#ffffff', '#eeaa77',usedoc)
        file.write(page(self.name,contents))
        file.close()

comment = None
param = None
params = None
stringpar = None
synopsis = None

def link(name):
    # for improved readability, 'sf' program name prefix will not be blue and underlined
    if name[0:2] == 'sf':
        prefix = '<font color="#778899">sf</font>'
        inlink = name[2:]
    else:
        prefix = ''
        inlink = name
    return '%s<a href="%s.html" title="%s">%s</a>' %\
    (prefix, progs[name].name, progs[name].desc, inlink)

def html(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    file = open (os.path.join(dir,'index.html'),'w')
    name = '<big><big><strong>Madagascar Programs</strong></big></big>'
    content = heading(name,'#ffffff','#7799ee')
    # Read subversion version number, if possible
    proc = os.popen('svn stat -v SConstruct 2>/dev/null')
    raw  = proc.read()
    proc.close()
    if len(raw) > 0: # SConstruct is under version control
        try:
            revnr = raw[13:].rsplit()[0]
            know_revnr = True
        except: # Python < 2.4 or strange svn output
            know_revnr = False
    else:
        know_revnr = False
    if have_datetime_module or know_revnr:
        content += 'Generated'
        if have_datetime_module:
            content += ' on ' + str(datetime.date.today())
        if know_revnr:
            content += ' from development version r' + revnr
    dirs = {}
    for prog in progs.keys():
        dir = os.path.dirname(progs[prog].file)
        if not dirs.has_key(dir):
            dirs[dir] = []
        dirs[dir].append(prog)
    keys = dirs.keys()
    keys.sort()
    for dir in keys:
        names = dirs[dir]
        names.sort()
        content = content + bigsection(dir,'#ffffff', '#ee77aa',
                                       multicolumn(names,link))
    file.write(page('all programs',content))
    file.close()

def text(dir,name):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    file = open (os.path.join(dir,name),'w')
    file.write('Madagascar Programs\n')
    dirs = {}
    for prog in progs.keys():
        dir = os.path.dirname(progs[prog].file)
        if not dirs.has_key(dir):
            dirs[dir] = []
        dirs[dir].append(prog)
    keys = dirs.keys()
    keys.sort()
    for dir in keys:
        names = dirs[dir]
        names.sort()
        file.write('\n[%s]\n\n%s\n' % (dir,'\n'.join(names)))
    file.close()

def spec(dir,name='extend.spec'):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    file = open (os.path.join(dir,name),'w')
    file.write('''
# Madagascar Enumerations

[enum-bool]
Desc:	boolean values
y	yes
n	no

[enum-file]
#Desc:	Known file/data types
rsf	Madagascar Regularly Sampled Format file
vpl	(unknown description)

[enum-color]
#Desc:	Port colors (#RRGGBB) to assign to file types
rsf	#5f5fc0
vpl	#c05f5f
''')

    # TO DO: Define tksu sections from most common parameters
    file.close()


# regular expressions
comment = {}
param = {}
params = {}
param2 = {}
params2 = {}
stringpar = {}
synopsis = {}
inpout = {}
version = {}

comment['python'] = re.compile(r'[^\'\"]*[\'\"]+([^\'\"]+)')
param['python'] = re.compile(r'par\.(?P<type>bool|int|float|string)'
                             '\s*\(\s*[\"\'](?P<name>\w+)[\"\']\s*'
                             '(?:\,\s*(?P<default>[^\)]+))?\)'
                             '(?:\s*\#\s*(?P<range>[\[][^\]]+[\]])?\s*'
                             '(?P<desc>[^#\n]+\S))?')
synopsis['python'] = re.compile(r'\s*\#\s*Takes\s*\:\s*'
                                '((?:[^\n]|\n\#[^\n])+)((?:.|\n)*)$')
inpout['python'] = re.compile(r'\s*(?P<name>\w+)\s*=\s*'
                                'rsf\.(?P<io>Input|Output)'
                                '\s*\(\s*(?:[\"\'](?P<tag>\w+)[\"\'])?')
version['python'] = re.compile(r'\#\s*\$Id\:\s*(.+\S)\s*\$/')

comment['f90'] = re.compile(r'(?:\!([^!\n]+)\n)+')
param['f90'] = re.compile(r'from_par\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                          '\s*(?P<var>[\w\_]+)\s*'
                          '(?:\,\s*(?P<default>[^\)]+))?\)'
                          '(?:\s*\!\s*(?P<range>[\[][^\]]+[\]])?\s*'
                          '(?P<desc>[^!\n]+\S))?')
synopsis['f90'] = re.compile(r'\s*\!\s*Takes\s*\:\s*'
                             '((?:[^\n]|\n\![^\n])+)((?:.|\n)*)$')
inpout['f90'] = re.compile(r'\s*(?P<name>\w+)\s*=\s*'
                    'rsf_(?P<io>input|output)'
                    '\s*\(\s*(?:\"(?P<tag>\w+)\")?')
version['f90'] = re.compile(r'\!\s*\$Id\:\s*(.+\S)\s*\$/')

comment['c'] = re.compile(r'\/\*((?:[^*]+|\*[^/])+)\*\/')
param['c'] = re.compile(r'(?:if\s*\(\!)?\s*sf_get'
                        '(?P<type>bool|largeint|int|float|double)'
                        '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                        '\s*\&(?P<var>[\w\_\[\]]+)\s*[\)]\s*[\)]?\s*'
                        '(?:[\{]|' # either \{ or
                        '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|'
                        'sf_[^\;]+)?' # or sf_error
                        '[\;])\s*' # ending with ;
                        '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                        '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
params['c'] = re.compile(r'sf_get(?P<type>bools|ints|floats|strings)'
                         '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                         '\s*(?P<var>[\w\_\[\]]+)\s*\,'
                         '\s*(?P<size>[\w\_]+)\s*\)\s*'
                         '[^\;\{]*[\;\{]\s*' # ending with ; or {
                         '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                         '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
param2['c'] = re.compile(r'sf_get(?P<type>bool|largeint|int|float|string)\s*'
                    '\([^/]+\/\*\(\s*(?P<name>[\w\#]+)'
                    '(?:=(?P<default>\S+))?'
                    '\s*(?P<desc>[^\)]+)\)\*\/')
params2['c'] = re.compile(r'sf_get(?P<type>bools|ints|floats|strings)'
                     '\s*\([^\,]+\,[^\,]+\,'
                     '\s*(?P<size>[\w\_]+)\s*\)[^/]+'
                     '\/\*\(\s*(?P<name>[\w\#]+)'
                     '(?:=(?P<default>\S+))?'
                     '\s*(?P<desc>(?:[^\)]|\)[^\*])+)\)\*\/')
stringpar['c'] = re.compile(r'sf_getstring\s*\(\s*\"(?P<name>\w+)\"'
                       '[^\;\{]*[\;\{]\s*(?:\/\*'
                       '\s*(?P<desc>(?:[^*]|\*[^/])+)\*\/)?')
synopsis['c'] = re.compile(r'\s*Takes\s*\:\s*((?:[^\n]|[\n][^\n])+)'
                      '((?:.|\n)*)$')
inpout['c'] = re.compile(r'\s*(?P<name>\w+)\s*=\s*'
                    'sf_(?P<io>input|output)'
                    '\s*\(\s*\"(?P<tag>\w+)\"')
version['c'] = re.compile(r'\/\*\s*\$Id\:\s*(.+\S)\s*\$\s*\*\/')


def getprog(file,out,lang = 'c',rsfprefix = 'sf',rsfsuffix='rsf',
            rsfplotprefix='vp',rsfplotsuffix='vpl'):
    global comment, param, params, param2, params2, \
           synopsis, stringpar, inpout, version
    name = rsfprefix + re.sub('^[MX]','',os.path.basename(file))
    if lang[:2] == 'py':
        name = re.sub('\.py$','',name)
    elif lang[0] =='f':
        name = re.sub('\.f\d*$','',name)
    else:
        name = re.sub('\.c(c|u)?$','',name)
    src = open(file,"r")   # open source
    text = ''.join(src.readlines())
    src.close()
    first = comment[lang].match(text)
    if first:
        tops = first.group(1).split("\n")
        desc = tops.pop(0).lstrip()
        first = '\n'.join(tops)
    else:
        desc = None
    prog = rsfprog(name,file,desc)
    file = re.sub('^[^\/]*\/','',file)
    out.write("%s = rsfdoc.rsfprog('%s','%s','''%s''')\n" %
              (name,name,file,desc))
    files = inpout[lang].findall(text)
    snps = name
    valid = {}
    for par in files:
        filename = par[0]
        if valid.get(filename,1):
            io = par[1]
            tag = par[2]
            if tag == 'in' or (not tag and io == 'input'):
                iostring = ' < %s.%s' % (filename,rsfsuffix)
            elif tag == 'out' or (not tag and io == 'output'):
                iostring = ' > %s.%s' % (filename,rsfsuffix)
            else:
                iostring = ' %s=%s.%s' % (tag,filename,rsfsuffix)
                type = 'file   '
                desc = 'auxiliary %s file name' % io
                prog.par(tag,rsfpar(type,desc=desc))
                out.write("%s.par('%s',rsfdoc.rsfpar('%s',desc='''%s'''))\n" %
                          (name,tag,type,desc))
            snps = snps + iostring
        valid[filename]=0
    if re.match(rsfplotprefix,name):
        snps = snps + ' > plot.' + rsfplotsuffix
    parline = ''
    pars = params.get(lang)
    if pars:
        for par in pars.findall(text):
            type = par[0]
            parname = par[1]
            size = par[3]
            range = par[4]
            desc = par[5] + ' [%s]' % size

            prog.par(parname,rsfpar(type,None,range,desc))
            out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                      (name,parname,type,'',range,desc))
            parline = parline + " %s=" % parname
    pars = param[lang].findall(text)
    for par in pars:
        if lang == 'python':
            type = par[0]
            parname = par[1]
            default = par[2]
            range = par[3]
            desc = par[4]
        elif lang == 'f90':
            type = ''
            parname = par[0]
            default = par[2]
            range = par[3]
            desc = par[4]
        else: # c
            type = par[0]
            parname = par[1]
            default = par[3]
            range = par[4]
            desc = par[5]

        if type == 'bool':
            if default == 'true' or default == 'True':
                default = 'y'
            elif default == 'false' or default == 'False':
                default = 'n'

        prog.par(parname,rsfpar(type,default,range,desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                  (name,parname,type,default,range,desc))
        parline = parline + " %s=%s" % (parname,default)
    pars = params2.get(lang)
    if pars:
        for par in pars.findall(text):
            type = par[0]
            size = par[1]
            parname = par[2]
            default = par[3]
            range = '' # for now
            desc = par[4] + ' [%s]' % size

            prog.par(parname,rsfpar(type,default,range,desc))
            out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                      (name,parname,type,default,range,desc))
            parline = parline + " %s=%s" % (parname,default)
    pars = param2.get(lang)
    if pars:
        for par in pars.findall(text):
            type = par[0]
            parname = par[1]
            default = par[2]
            range = '' # for now
            desc = par[3]

            prog.par(parname,rsfpar(type,default,range,desc))
            out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                      (name,parname,type,default,range,desc))
            parline = parline + " %s=%s" % (parname,default)
    pars = stringpar.get(lang)
    if pars:
        for par in pars.findall(text):
            type = 'string '
            parname = par[0]
            desc = par[1]
            if parname in prog.pars: # defined previously as file
                filedesc = prog.pars[parname].desc.strip()
                if desc:
                    desc = desc + '(%s)' % filedesc
                else:
                    desc = filedesc
            else:
                parline = parline + " %s=" % (parname)
            prog.par(parname,rsfpar(type,desc=desc))
            out.write("%s.par('%s',rsfdoc.rsfpar('%s',desc='''%s'''))\n" %
                      (name,parname,type,desc))
    snps = snps + parline
    base = name[len(rsfprefix):]
    if base in docprogs:
        wiki = r'http://reproducibility.org/wiki/Guide_to_madagascar_programs#sf'+base
        prog.weblink(wiki)
        out.write("%s.weblink('%s')\n" % (name,wiki))
    vers = version[lang].search(text)
    if vers:
        prog.version(vers.group(1))
        out.write("%s.version('''%s''')\n" % (name,vers.group(1)))
    if not first:
        first = ''
    else:
        info = synopsis[lang].match(first)
        if info:
            snps = snps + ' ' + info.group(1).lstrip()
            first = info.group(2).lstrip()
    prog.synopsis(snps,first)
    out.write("%s.synopsis('''%s''','''%s''')\n" % (name,snps,first))
    out.write("rsfdoc.progs['%s']=%s\n\n" % (name,name))

def cli(rsfprefix = 'sf',rsfplotprefix='vp'):
    # Implements own UI instead of Madagascar's standard Python API
    # for UI compatibility with pydoc
    import getopt
    import rsfprog

    this = sys.argv.pop(0)
    class BadUsage: pass

    try:
        opts, args = getopt.getopt(sys.argv, 'k:w:t:s:m:g:l:r:')
        dir = None
        typ = None
        rep = os.getcwd()
        for opt, val in opts:
            if opt[0] == '-' and opt[1] in 'wtsmgl':
                dir = val
                typ = opt[1]
            elif opt == '-r':
                rep = val
            elif opt == '-k':
                val = val.lower()
                doc = ''
                for prog in progs.keys():
                    desc = progs[prog].desc
                    if re.search(val,desc.lower()):
                        doc = doc + "%s: %s\n" % (bold(prog),desc)
                pydoc.pager(doc)
                return

        if not args:
            if typ == 'w':
                html(dir)
                for prog in progs.keys():
                    main = progs.get(prog)
                    if main:
                        main.html(dir,rep)
            elif typ == 't':
                text(dir,'INDEX.txt')
                for prog in progs.keys():
                    main = progs.get(prog)
                    if main:
                        main.text(dir)
            elif typ == 's':
                spec(dir,'extend.spec')
                for prog in progs.keys():
                    main = progs.get(prog)
                    if main:
                        main.spec(dir)
            elif typ == 'g':
                text(dir,'index.man')
                for prog in progs.keys():
                    main = progs.get(prog)
                    if main:
                        main.man(dir)
            else:
                raise BadUsage

        for prog in args:
            if prog == 'vppen' or \
                    not re.match(rsfprefix,prog) and \
                    not re.match(rsfplotprefix,prog):
                prog = rsfprefix + prog
            main = progs.get(prog)
            if main:
                if   typ == 'w':
                    main.html(dir,rep)
                elif typ == 't':
                    main.text(dir,prog)
                elif typ == 's':
                    main.spec(dir,prog)
                elif typ == 'm':
                    main.mwiki(dir,prog)
                elif typ == 'g':
                    main.man(dir,prog)
                elif typ == 'l':
                    main.latex(dir,prog)
                else:
                    main.document()
            else:
                print "No program %s in Madagascar." % prog

    except (getopt.error, BadUsage):
        print '''sfdoc - the RSF documentation tool

%(prog)s <prog1> <prog2> ...
    Show documentation on programs.

%(prog)s -t <dir> [<prog1> <prog2> ... ]
    Write plain text documentation in <dir> directory

%(prog)s -s <dir> [<prog1> <prog2> ... ]
    Write TKSU block specs in <dir> directory

%(prog)s -w <dir> [-r <rep>] [<prog1> <prog2> ...]
    Write program HTML documentaton in <dir> directory, optional <rep> referes to repository.

%(prog)s -m <dir>  <prog1> <prog2> ...
    Write program documentaton in MediaWiki format in <dir> directory.

%(prog)s -g <dir>  <prog1> <prog2> ...
    Write program documentaton in groff (man page) format in <dir> directory.

%(prog)s -l <dir> <prog1> <prog2> ...
    Write program LaTeX documentaton in <dir> directory.

%(prog)s -k <keyword>
    Search for a keyword in the description lines of all available programs.
''' % {'prog':this}

if __name__ == "__main__":
    junk = open('junk.py',"w")
    junk.write("import rsfdoc\n\n")
    junk.write("rsfprog = {}\n")
    getprog('filt/main/dd.c',junk)
    junk.write("sfdd.document()\n\n")
    junk.close()
    #
    import junk
    #
    os.unlink("junk.py")
    os.unlink("junk.pyc")
