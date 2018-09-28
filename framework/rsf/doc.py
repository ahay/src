# Copyright (C) 2004 University of Texas at Austin
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




# this produces errors with kdmig2d.py - the first su program doc
#from __future__ import unicode_literals

from __future__ import division, absolute_import, print_function
import pydoc, re, sys, os, glob, signal
import rsf.path as rsfpath

try:
    import datetime
    have_datetime_module = True
except: # Python < 2.3
    have_datetime_module = False

progs = {}
data = {}

# Programs with wiki documentation -- this should be parsed from a downloaded
# m8r wiki page, with fallback info in a config file

docprogs = '''
add attr awefd2d cat cmplx conjgrad cp csv2rsf cut dd disfil dottest get
headerattr headercut headermath headersort headerwindow in interleave mask math
pad pick prep4plot put real remap1 reverse rm rotate rtoc scale segyread
segywrite spike spray srmig3 stack stretch transp window sizes figlist booklist
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
    child = os.fork()
    if child:
        (pid,exit) = os.waitpid(child,0)
        child = 0
        return exit
    else:
        os.system(comm)
        os._exit(0)

def subdirs():
    return [x for x in list(filter(os.path.isdir,glob.glob('[a-z]*'))) if x[-5:] != '_html']

def getversion(target=None,source=None,env=None):
    label = env.get('version')
    name = str(target[0])
    out = open(name,'w')
    out.write('#!/usr/bin/env python\n\n')
    out.write('label="%s"\n\n' % label)
    out.write('''
if __name__ == "__main__":
    print(label)
    ''')
    out.close()
    os.chmod(name,0o775)
    return 0

def getprogs(target=None,source=None,env=None):
    out = open(str(target[0]),'w')
    out.write('import sys, os\n\n')
    out.write('import rsf.doc\n\n')
    for mod in map(str,source):
        py = os.path.splitext(os.path.basename(mod))[0]
        out.write('import rsf.%s\n' % py)
    out.write('\nimport rsf.vpplot\n\n')
    out.write('''
try:
    import rsf.use
except:
    pass

RSFROOT="%s"

def selfdoc():
    'Display man page'
    prognm = os.path.basename(sys.argv[0])
    if prognm[0] == 'M' and prognm[-3:] == '.py':
        # User testing code in his local directory
        prognm = 'sf' + prognm[1:-3]
        msg = 'Self-doc may be out of synch, do "scons install" in RSFSRC'
        sys.stderr.write(msg+'\\n')

    prog = rsf.doc.progs.get(prognm)
    if prog != None:
        prog.document(25,RSFROOT)
    else:
        sys.stderr.write('No installed man page for ' + prognm+'\\n')
''' % env.get('RSFROOT'))
    out.close()

def use(target=None,source=None,env=None):
    'Collect program uses information'
    out = open(str(target[0]),'w')
    out.write('import rsf.doc\n\n')
    for dotproj in map(str,source):

        try:
            glo = {}
            loc = {}
            exec(compile(open(dotproj).read(), dotproj, 'exec'),glo,loc)
        except:
            sys.stderr.write('problem with %s' % dotproj)
            continue

        dirname = os.path.dirname(dotproj)
        dirname,project = os.path.split(dirname)
        dirname,chapter = os.path.split(dirname)
        dirname,book    = os.path.split(dirname)

        for prog in loc['uses']:
            out.write('rsf.doc.progs["%s"].use("%s","%s","%s")\n' %
                      (prog,book,chapter,project))
        out.write('\n')
    out.close()
    return 0

def selfdoc(target=None,source=None,env=None):
    src = str(source[0])
    doc = open(str(target[0]),"w")
    rsfprefix = env.get('rsfprefix','sf')
    rsfsuffix = env.get('rsfsuffix','rsf')
    lang = env.get('lang','c')
    known_version = env.get('version','unknown')
    strip = env.get('strip',None)
    getprog(src,doc,lang,rsfprefix,rsfsuffix,known_version=known_version,strip=strip)
    doc.close()
    return 0

def bold(text):
    """Format a string in bold by overstriking."""
    return ''.join([ch + "\b" + ch for ch in text])

def underline(text):
    """Format a string in underline by overstriking."""
    return ''.join([ch + "\b_" for ch in text])

def underline_match(text):
    """Underline a matching expression"""
    return underline(text.group(0))

def replace(text, *pairs):
    """Do a series of global replacements on a string."""
    while pairs:
        text = pairs[1].join(text.split(pairs[0]))
        pairs = pairs[2:]
    return text

def section(head,body):
    text = "\n".join(["\t" + line for line in body.split("\n")])
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
    return html_section(*(title,)+args)

def multicolumn(list, format, cols=4):
    """Format a list of items into a multi-column list."""
    result = ''
    rows = (len(list)+cols-1)//cols
    for col in range(cols):
        result = result + '<td width="%d%%" valign=top>' % (100//cols)
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
            elif (desc.find('Output')>=0):
                rw = 'w'
                if not name:
                    name = 'out'
            elif (desc.find('Input')>=0):
                rw = 'r'
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
            (self.type,name,self.default,self.range,self.desc.rstrip())
        return troff
class rsfdata(object):
    def __init__(self,name):
        self.name = name
        self.uses = {}
    def use (self,book,chapter,project):
        if book not in self.uses:
            self.uses[book]={}
        if chapter not in self.uses[book]:
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
        if book not in self.uses:
            self.uses[book]={}
        if chapter not in self.uses[book]:
            self.uses[book][chapter] = []
        self.uses[book][chapter].append(project)
    def docstring(self,usedoc_max,rsfroot):
        doc = section('name',self.name)
        if self.desc:
            doc = doc + section('description',self.desc)
        if self.snps:
            doc = doc + section('synopsis',self.snps)
        if self.cmts:
            cmts = re.sub(r'http://[\S]+',underline_match,self.cmts)
            doc = doc + section('comments',cmts)
        pars =  list(self.pars.keys())
        if pars:
            pars.sort()
            pardoc = ''
            for par in pars:
                pardoc = pardoc + self.pars[par].show(par)
            doc = doc + section('parameters',pardoc.rstrip())
        if self.also:
            doc = doc + section('see also',self.also)
        books = list(self.uses.keys())
        if books:
            usedoc = ''
            usedoc_i = 0
            books.sort()
            for book in books:
                chapters = list(self.uses[book].keys())
                chapters.sort()
                for chapter in chapters:
                    for project in self.uses[book][chapter]:
                        if usedoc_i < usedoc_max:
                            usedoc += '%s/%s/%s\n' % (book,chapter,project)
                        usedoc_i += 1
            if usedoc:
                if usedoc_i > usedoc_max:
                    usedoc += '%d more examples listed in:\n' % \
                              (usedoc_i - usedoc_max)
                    usedoc += '%s/share/doc/madagascar/html/%s.html\n'%\
                              (rsfroot,self.name)
            doc = doc + section('used in',usedoc.rstrip())
        doc = doc + section('source',self.file)
        if self.wiki:
            doc = doc + section('documentation',underline(self.wiki))
        if self.vers:
            doc = doc + section('version',self.vers)
        return doc
    def document(self,usedoc_max,rsfroot):
        doc = self.docstring(usedoc_max,rsfroot)
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
        pars =  list(self.pars.keys())
        if pars:
            pars.sort()
            for par in pars:
                contents = contents + self.pars[par].mwiki(par)
        contents = contents+'|}\n'
        file.write(contents)
        file.close()
    def man(self,dir,usedoc_max,rsfroot,name=None):
        if not name:
            name = self.name
        file = open (os.path.join(dir,name + '.1'),'w', encoding='utf-8')
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
        pars =  list(self.pars.keys())
        if pars:
            contents += '.SH PARAMETERS\n.PD 0\n'
            pars.sort()
            for par in pars:
                contents += self.pars[par].man(par)
        books = list(self.uses.keys())
        if books:
            usedoc = ''
            usedoc_i = 0
            books.sort()
            for book in books:
                chapters = list(self.uses[book].keys())
                chapters.sort()
                for chapter in chapters:
                    for project in self.uses[book][chapter]:
                        if usedoc_i < usedoc_max:
                            usedoc += '.TP\n.I %s/%s/%s\n' % \
                                      (book,chapter,project)
                        usedoc_i += 1
            if usedoc:
                contents = contents + '.SH USED IN\n%s' % usedoc
                if usedoc_i > usedoc_max:
                    contents += '.TP\n%d more examples listed in:\n' % \
                                (usedoc_i - usedoc_max)
                    contents += '.TP\n%s/share/doc/madagascar/html/%s.html\n'%\
                                 (rsfroot,name)
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
        pars =  list(self.pars.keys())
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
        if dir:
            file = open (os.path.join(dir,name + '.txt'),'w', encoding='utf-8')
        contents = 'Program %s | %s\n' % (name,self.desc)
        if self.snps:
            contents = contents + '[SYNOPSIS]\n%s\n' % self.snps
        if self.cmts:
            contents = contents + '[COMMENTS]\n%s\n' % self.cmts
        pars =  list(self.pars.keys())
        if pars:
            contents = contents + '[PARAMETERS]\n'
            #sys.stderr.write('type(pars)=%s\n'%type(pars))
            #sys.stderr.write('pars=%s\n'%repr(pars))
            #pars.sort()
            for par in sorted(pars):
                contents = contents + self.pars[par].text(par)
        filedir = os.path.split(self.file)[0]
        if filedir:
            contents = contents + '[DIRECTORY]\n%s\n' % filedir
        if dir:
            file.write(contents)
            file.close()
        return contents
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
                    data, ext  = filename.split('.')
                except:
                    continue
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
        pars =  list(self.pars.keys())
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
        hfile = open (os.path.join(dir,self.name + '.html'),'w')
        name = '<big><big><strong>%s</strong></big></big>' % self.name
        if self.vers:
            name = name + " (%s)" % self.vers
        if self.wiki:
            wiki = '<br><a href="%s">Documentation</a>' % self.wiki
        else:
            wiki = ''
        contents = heading(name,'#ffffff','#7799ee',
                           '<a href="./index.html">index</a><br>'+
                           '<a href="%s/%s">%s</a>%s' %
                           (rep,self.file,self.file,wiki))
        if self.desc:
            contents += self.desc
        if self.snps:
            contents += bigsection('Synopsis','#fffff', '#aa55cc',self.snps)
        if self.cmts:
            cmts = re.sub(r'(http://[\S]+)',r'<a href="\1">\1</a>',self.cmts)
            contents += cmts.replace('\n','<br>\n')
        pars =  list(self.pars.keys())
        if pars:
            pars.sort()
            pardoc = ''
            bgcol = '#ffc8d8'
            for par in pars:
                pardoc += heading(self.pars[par].html(par),'#000000', bgcol,
                                  self.pars[par].desc.replace('\n','<br>\n'),
                                  add='') + '\n'
                if bgcol=='#ffc8d8':
                    bgcol ='#f0f0f8'
                else:
                    bgcol = '#ffc8d8'
            contents += bigsection('Parameters','#ffffff', '#ee77aa',pardoc)
        books = list(self.uses.keys())
        if books:
            usedoc = ''
            books.sort()
            for book in books:
                bookdoc = ''
                chapters = list(self.uses[book].keys())
                chapters.sort()
                for chapter in chapters:
                    for project in self.uses[book][chapter]:
                        proj = os.path.join(chapter,project)
                        bookdoc += '''
                        <a href="book/%s/%s.html">%s</a><br>
                        ''' % (book,proj,proj)
                usedoc += bigsection(book.upper(),'#000000','#ffd8c8',bookdoc)
            contents += bigsection('Used In','#ffffff', '#eeaa77',usedoc)
        hfile.write(page(self.name,contents.encode('utf-8')))
        hfile.close()

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

def html(dir,known_version):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    file = open (os.path.join(dir,'index.html'),'w')
    name = '<big><big><strong>Madagascar Programs</strong></big></big>'
    content = heading(name,'#ffffff','#7799ee')

    if known_version[-4:] == '-svn' or known_version == '':
        # Read subversion version number, if possible
        proc = os.popen('svn stat -v SConstruct 2>/dev/null')
        raw  = proc.read()
        proc.close()
        if len(raw) > 0: # file is under version control
            try:
                revnr = raw[13:].rsplit()[0]
                known_version += ' r'+revnr
            except: # Python < 2.4 or strange svn output
                pass

    rev_add2content = 'version ' + known_version

    if have_datetime_module or know_revnr:
        content += 'Generated'
        if have_datetime_module:
            content += ' on ' + str(datetime.date.today())
        if rev_add2content != '':
            content += ' from ' + rev_add2content
    dirs = {}
    for prog in list(progs.keys()):
        dir = os.path.dirname(progs[prog].file)
        if dir not in dirs:
            dirs[dir] = []
        dirs[dir].append(prog)
    keys = list(dirs.keys())
    keys.sort()
    for dir in keys:
        names = dirs[dir]
        names.sort()
        content = content + bigsection(dir,'#ffffff', '#ee77aa',
                                       multicolumn(names,link))
    file.write(page('all programs',content.encode('utf-8')))
    file.close()

def text(dir,name):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    file = open (os.path.join(dir,name),'w')
    file.write('Madagascar Programs\n')
    dirs = {}
    for prog in list(progs.keys()):
        dir = os.path.dirname(progs[prog].file)
        if dir not in dirs:
            dirs[dir] = []
        dirs[dir].append(prog)
    keys = list(dirs.keys())
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
chars = {}

comment['python'] = re.compile(r'[^\'\"]*([\'\"]+)(?P<comment>[^\'\"].+?)\1',
                               re.DOTALL)
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
chars['python'] = r'# '

comment['c++'] = re.compile(r'\/\/(?P<comment>[^\n]+)\n')
inpout['c++'] = re.compile(r'(?P<io>iRSF|oRSF)\s*(?P<name>\w+)'
                           '\s*(?:\(\s*\"?(?P<tag>\w+)\"?\s*\))?')
param['c++'] = re.compile(r'par.get\s*\(\s*[\"\'](?P<name>\w+)[\"\']\s*'
                          '(?:\,\s*(?P<var>[^\)\,]+))'
                          '(?:\,\s*(?P<default>[^\)]+))?\)'
                          '(?:\s*\;\s*\/\/\s*(?P<range>[\[][^\]]+[\]])?\s*'
                          '(?P<desc>[^\n]+\S))?')
version['c++'] = re.compile(r'\/\/\s*\$Id\:\s*(.+\S)\s*\$/')
chars['c++'] = r'/ '

comment['f90'] = re.compile(r'\!(?P<comment>[^\n]+(?:\n\![^\n]+)+)\n')
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
chars['f90'] = r'! '

comment['c'] = re.compile(r'\/\*(?P<comment>(?:[^*]+|\*[^/])+)\*\/')
param['c'] = re.compile(r'(?:if\s*\(\!)?\s*sf_get'
                        '(?P<type>bool|largeint|int|float|double)'
                        '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                        '\s*\&(?P<var>[\w\_\[\]\.\+\-\>]+)\s*[\)]\s*[\)]?\s*'
                        '(?:[\{]|' # either \{ or
                        '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|'
                        'sf_[^\;]+)?' # or sf_error
                        '[\;])\s*' # ending with ;
                        '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                        '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
params['c'] = re.compile(r'sf_get(?P<type>bools|ints|floats|strings)'
                         '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                         '\s*(?P<var>[\w\_\[\]\.]+)\s*\,'
                         '\s*(?P<size>[\w\_\-\+]+)\s*\)\s*'
                         '[^\;\{]*[\;\{]\s*' # ending with ; or {
                         '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                         '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
param2['c'] = re.compile(r'sf_get(?P<type>bool|largeint|int|float|string)\s*'
                    '\((?:[^/]|/[^\*])+\/\*\(\s*(?P<name>[\w\#]+)'
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
chars['c'] = ' '


def getprog(file,out,lang = 'c',rsfprefix = 'sf',rsfsuffix='rsf',
            rsfplotprefix='vp',rsfplotsuffix='vpl',known_version='unknown',strip = None):
    global comment, param, params, param2, params2, \
           synopsis, stringpar, inpout, version
    name = re.sub('^[MX_]','',os.path.basename(file))
    if strip: # remove leading characters in the file name
        source = os.path.join(os.path.dirname(file),name)
    else:
        source = file
    name = rsfprefix + name
    if lang[:2] == 'py':
        name = re.sub('\.py$','',name)
    elif lang[0] =='f':
        name = re.sub('\.f\d*$','',name)
    else:
        name = re.sub('\.c(c|u)?$','',name)
    cname = re.sub('\-','',name)
    src = open(file,"r")   # open source
    text = ''.join(src.readlines())
    src.close()
    first = comment[lang].match(text)
    if first:
        tops = first.group('comment').split("\n")
        desc = tops.pop(0).lstrip()
        first = '\n'.join([x.lstrip(chars[lang]) for x in tops])
    else:
        desc = None
    prog = rsfprog(name,source,desc)
    source = re.sub('^[^\/]*\/','',source)
    out.write("%s = rsf.doc.rsfprog('%s','%s','''%s''')\n" %
              (cname,name,source,desc))
    files = inpout[lang].findall(text)
    snps = name
    valid = {}
    for par in files:
        if lang == 'c++':
            io = par[0]
            filename = par[1]
            tag = par[2]
            if tag == '0':
                valid[filename]=0
        else:
            filename = par[0]
            io = par[1]
            tag = par[2]
        if valid.get(filename,1):
            iotype = {'i':'input','o':'output'}[io[0].lower()]
            if tag == 'in' or (not tag and iotype == 'input'):
                iostring = ' < %s.%s' % (filename,rsfsuffix)
            elif tag == 'out' or (not tag and iotype == 'output'):
                iostring = ' > %s.%s' % (filename,rsfsuffix)
            else:
                iostring = ' %s=%s.%s' % (tag,filename,rsfsuffix)
                type = 'file   '
                desc = 'auxiliary %s file name' % iotype
                prog.par(tag,rsfpar(type,desc=desc))
                out.write("%s.par('%s',rsf.doc.rsfpar('%s',desc='''%s'''))\n" %
                          (cname,tag,type,desc))
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
            out.write("%s.par('%s',rsf.doc.rsfpar('%s','%s','%s','''%s'''))\n" %
                      (cname,parname,type,'',range,desc))
            parline = parline + " %s=" % parname
    pars = param[lang].findall(text)
    for par in pars:
        if lang == 'python':
            type = par[0]
            parname = par[1]
            default = par[2]
            if type == 'string':
                default = default.strip('"').strip("'")
            range = par[3]
            desc = par[4]
        elif lang == 'f90':
            type = ''
            parname = par[0]
            default = re.sub('.true.','y',re.sub('.false.','n',par[2]))
            range = par[3]
            desc = par[4]
        elif lang == 'c++':
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
            bool_default = default.lower().strip()
            if bool_default == 'true':
                default = 'y'
            elif bool_default == 'false':
                default = 'n'

        prog.par(parname,rsfpar(type,default,range,desc))
        out.write("%s.par('%s',rsf.doc.rsfpar('%s','%s','%s','''%s'''))\n" %
                  (cname,parname,type,default,range,desc))
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
            out.write("%s.par('%s',rsf.doc.rsfpar('%s','%s','%s','''%s'''))\n" %
                      (cname,parname,type,default,range,desc))
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
            out.write("%s.par('%s',rsf.doc.rsfpar('%s','%s','%s','''%s'''))\n" %
                      (cname,parname,type,default,range,desc))
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
            out.write("%s.par('%s',rsf.doc.rsfpar('%s',desc='''%s'''))\n" %
                      (cname,parname,type,desc))
    snps = snps + parline
    base = name[len(rsfprefix):]
    if base in docprogs:
        wiki = r'http://ahay.org/wiki/Guide_to_madagascar_programs#sf'+base
        prog.weblink(wiki)
        out.write("%s.weblink('%s')\n" % (cname,wiki))
    vers = version[lang].search(text)
    if vers:
        known_version += ' ' + vers.group(1)
    prog.version(known_version)
    out.write("%s.version('%s')\n" % (cname, known_version))
    if not first:
        first = ''
    else:
        info = synopsis[lang].match(first)
        if info:
            snps = snps + ' ' + info.group(1).lstrip()
            first = info.group(2).lstrip()
    prog.synopsis(snps,first)
    out.write("%s.synopsis('''%s''','''%s''')\n" % (cname,snps,first))
    out.write("rsf.doc.progs['%s']=%s\n\n" % (name,cname))

def cli(rsfprefix = 'sf',rsfplotprefix='vp'):
    # Implements own UI instead of Madagascar's standard Python API
    # for UI compatibility with pydoc
    import getopt
    import rsf.prog

    this = sys.argv.pop(0)
    class BadUsage: pass

    root = rsf.prog.RSFROOT

    try:
        opts, args = getopt.getopt(sys.argv, 'k:w:t:s:m:g:l:r:v:u:')
        dir = None
        typ = None
        rep = os.getcwd()
        ver = 'unknown'
        usedoc_max = 25
        for opt, val in opts:
            if opt[0] == '-' and opt[1] in 'wtsmgl':
                dir = val
                typ = opt[1]
            elif opt == '-r':
                rep = val
            elif opt == '-v':
                ver = val
            elif opt == '-u':
                usedoc_max = int(val)
            elif opt == '-k':
                val = val.lower()
                doc = ''
                for prog in list(progs.keys()):
                    desc = progs[prog].desc
                    if re.search(val,desc.lower()):
                        doc = doc + "%s: %s\n" % (bold(prog),desc)
                pydoc.pager(doc)
                return

        if not args:
            if typ == 'w':
                html(dir,ver)
                for prog in list(progs.keys()):
                    main = progs.get(prog)
                    if main:
                        main.html(dir,rep)
            elif typ == 't':
                text(dir,'INDEX.txt')
                for prog in list(progs.keys()):
                    main = progs.get(prog)
                    if main:
                        main.text(dir)
            elif typ == 's':
                spec(dir,'extend.spec')
                for prog in list(progs.keys()):
                    main = progs.get(prog)
                    if main:
                        main.spec(dir)
            elif typ == 'g':
                text(dir,'index.man')
                for prog in list(progs.keys()):
                    main = progs.get(prog)
                    if main:
                        main.man(dir,usedoc_max,root)
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
                    main.man(dir,usedoc_max,root,prog)
                elif typ == 'l':
                    main.latex(dir,prog)
                else:
                    main.document(usedoc_max,root)
            else:
                print ('''Run %s with parameters.
To obtain a selfdoc, install %s with Madagascar: http://www.ahay.org/wiki/Adding_new_programs_to_Madagascar
                       ''' % (prog,prog))

    except (getopt.error, BadUsage):
        print('''sfdoc - the RSF documentation tool

%(prog)s [-u 25] <prog1> <prog2> ...
    Show documentation on programs.

%(prog)s -t <dir> [<prog1> <prog2> ... ]
    Write plain text documentation in <dir> directory

%(prog)s -s <dir> [<prog1> <prog2> ... ]
    Write TKSU block specs in <dir> directory

%(prog)s -w <dir> [-r <rep>] [-v <ver>] [<prog1> <prog2> ...]
    Write program HTML documentaton in <dir> directory, optional <rep> refers to repository, optional <ver> refers to version.

%(prog)s -m <dir> <prog1> <prog2> ...
    Write program documentaton in MediaWiki format in <dir> directory.

%(prog)s -g <dir> [-u 25]  <prog1> <prog2> ...
    Write program documentaton in groff (man page) format in <dir> directory.

%(prog)s -l <dir> <prog1> <prog2> ...
    Write program LaTeX documentaton in <dir> directory.

%(prog)s -k <keyword>
    Search for a keyword in the description lines of all available programs.
              ''' % {'prog':this})

if __name__ == "__main__":
    junk = open('junk.py',"w")
    junk.write("import rsf.doc\n\n")
    junk.write("rsfprog = {}\n")
    getprog('filt/main/dd.c',junk)
    junk.write("sfdd.document()\n\n")
    junk.close()
    #
    import junk
    #
    os.unlink("junk.py")
    os.unlink("junk.pyc")
