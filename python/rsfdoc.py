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
import re, sys, os, string, glob, string

progs = {}
data = {}

def subdirs():
    return filter(lambda x: x[-5:] != '_html',
                  filter(os.path.isdir,glob.glob('[a-z]*')))

def getprogs(target=None,source=None,env=None):
    out = open(str(target[0]),'w')
    dirs = env.get('dirs')
    out.write('import rsfdoc\n\n')
    for mod in dirs:
        out.write('import sf%s\n' % mod)
    out.write('\nimport vpplot\n\nimport rsfuse\n')
    out.close()

def use(target=None,source=None,env=None):
    out = open(str(target[0]),'w')
    doc = ['import rsfdoc\n']
    cwd = os.getcwd()
    bookdir = env.get('book')
    if os.path.isdir(bookdir):
        os.chdir(bookdir)
        for book in subdirs():
            os.chdir(book)
            print "%s..." % book
            for chapter in subdirs():
                os.chdir(chapter)
                print "...%s" % chapter
                for project in subdirs():
                    os.chdir(project)

                    #  (status,progs) = \
                    #  commands.getstatusoutput('scons -s .sf_uses')

                    for case in ('uses',):
                        sin, sout, serr = os.popen3('scons -s .sf_'+case)
                        sin.close()
                        progs = sout.read()
                        sout.close()
                        status = serr.read()
                        serr.close() 

                        if status:
                            print ('No %s found in book/%s/%s/%s/: %s' %
                                   (case,book,chapter,project,status))
                        elif string.find(progs,'scons') < 0:
                            for prog in string.split(progs):
                                doc.append(
                                    'rsfdoc.%s["%s"].use("%s","%s","%s")' %
                                    (('progs','data')[case=='data'],
                                     prog,book,chapter,project))
                    os.chdir('..')
                os.chdir('..')
            os.chdir('..')
        os.chdir(cwd)
    out.write(string.join(doc,'\n') + '\n')
    out.close()

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
    return string.join(map(lambda ch: ch + "\b" + ch, text),'')

def underline(text):
    """Format a string in underline by overstriking."""
    return string.join(map(lambda ch: ch + "\b_", text),'')

def replace(text, *pairs):
    """Do a series of global replacements on a string."""
    while pairs:
        text = string.join(string.split(text, pairs[0]), pairs[1])
        pairs = pairs[2:]
    return text

def section(head,body):
    text = string.join(map(lambda line: "\t" + line,
                           string.split(body,"\n")),"\n")
    return bold(string.upper(head)) + "\n" + text + "\n"

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

class rsfpar:
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
    def wiki(self,name):
        desc = string.replace(self.desc,'\n','<br>')
        desc = re.sub(r'<br>\s+(\w)',r'\n:\1',desc)
        desc = re.sub(r'(?:\s*<br>)+$','',desc)
        return "|-\n| ''%s'' || '''%s%s''' || %s || %s\n" % \
               (self.type,name,self.default,self.range,desc)
    def latex(self,name):
        tex = '\\underline{%s} & \\textbf{%s%s} & %s & ' % \
              (self.type,name,self.default,self.range)
        return tex + self.desc + '\\\\ \n'

class rsfdata:
    def __init__(self,name):
        self.name = name
        self.uses = {}
    def use (self,book,chapter,project):
        if not self.uses.has_key(book):
            self.uses[book]={}
        if not self.uses[book].has_key(chapter):
            self.uses[book][chapter] = []
        self.uses[book][chapter].append(project)
        
class rsfprog:
    def __init__(self,name,file,desc=None):
        self.name = name
        self.file = file
        self.desc = desc
        self.snps = None
        self.cmts = None
        self.also = None
        self.vers = None
        self.uses = {}
        self.pars = {}
    def check(self,par):
        '''par is name=value pair, check if it makes sense, return 1 if not'''
        pair = string.split(par,'=')
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
    def document(self):
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
            doc = doc + section('parameters',string.rstrip(pardoc))
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
            doc = doc + section('used in',string.rstrip(usedoc))
        doc = doc + section('source',self.file)
        if self.vers:
            doc = doc + section('version',self.vers)
        pydoc.pager(doc)
    def wiki(self,dir):
        file = open (os.path.join(dir,self.name + '.wiki'),'w')
        contents = '==%s==\n{| class="wikitable" align="center"\n' % \
                   self.name
        desc = '! colspan="4" style="background:#ffdead;" | %s\n' % self.desc
        contents = contents + desc
        if self.snps:
            contents = contents + '|-\n! colspan="4" | %s\n' % self.snps
        if self.cmts:
            cmts = string.replace(self.cmts,'\n','<br>')
            cmts = re.sub(r'(?:\s*<br>)+$','',cmts)
            contents = contents + '|-\n|  colspan="4" | %s\n' % cmts
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            for par in pars:
                contents = contents + self.pars[par].wiki(par)
        contents = contents+'|}\n'
        file.write(contents)
        file.close()
    def latex(self,dir):
        file = open (os.path.join(dir,self.name + '.tex'),'w')
        contents = '\\footnotesize\n'
        name = '\\subsection{%s: %s}\n' % (self.name,self.desc)
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
        contents = re.sub('_','\_',contents)
        file.write(contents)
        file.close()
    def html(self,dir,rep):
        file = open (os.path.join(dir,self.name + '.html'),'w')
        name = '<big><big><strong>%s</strong></big></big>' % self.name
        if self.vers:
            name = name + " (%s)" % self.vers
        contents = heading(name,'#ffffff','#7799ee',
                           '<a href="./index.html">index</a><br>'+
                           '<a href="%s/%s?view=markup">%s</a>' %
                           (rep,self.file,self.file))
        if self.desc:
            contents = contents + self.desc
        if self.snps:
            contents = contents + \
                       bigsection('Synopsis','#fffff', '#aa55cc',self.snps)
        if self.cmts:
            contents = contents + string.replace(self.cmts,'\n','<br>\n')
        pars =  self.pars.keys()
        if pars:
            pars.sort()
            pardoc = ''
            bgcol = '#ffc8d8'
            for par in pars:
                pardoc = pardoc + \
                         heading(self.pars[par].html(par),
                                 '#000000', bgcol,
                                 string.replace(self.pars[par].desc,
                                                '\n','<br>\n'),add='') + '\n'
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
                         bigsection(string.upper(book),
                                    '#000000','#ffd8c8',bookdoc)
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
    return '<a href="%s.html">%s</a>' % (progs[name].name, name)

def html(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)
    file = open (os.path.join(dir,'index.html'),'w')
    name = '<big><big><strong>RSF Programs</strong></big></big>'
    content = heading(name,'#ffffff','#7799ee')
    content = content + 'Add a search feature later.'
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
    file.write(page('RSF Programs',content))
    file.close()

def getprog(file,out,lang = 'c',rsfprefix = 'sf',rsfsuffix='rsf',
            rsfplotprefix='vp',rsfplotsuffix='vpl'):
    global comment, param, params, param2, params2, \
           synopsis, stringpar, inpout, version
    if not comment:
        if lang == 'f90':
            comment = re.compile(r'\!([^!]+)')
            param = re.compile(r'from_par\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                               '\s*(?P<var>[\w\_]+)\s*'
                               '(?:\,\s*(?P<default>[^\)]+))?\)' 
                               '(?:\s*\!\s*(?P<range>[\[][^\]]+[\]])?\s*'
                               '(?P<desc>[^!\n]+\S))?') # comment
            params = None
            param2 = None
            params2 = None
            stringpar = None
            synopsis = re.compile(r'\s*\!\s*Takes\s*\:\s*'
                                  '((?:[^\n]|\n\![^\n])+)((?:.|\n)*)$')
            inpout = re.compile(r'\s*(?P<name>\w+)\s*=\s*'
                                'rsf_(?P<io>input|output)'
                                '\s*\(\s*(?:\"(?P<tag>\w+)\")?')
            version = re.compile(r'\!\s*\$Id\:\s*(.+\S)\s*\$/')
        else: # c
            comment = re.compile(r'\/\*((?:[^*]+|\*[^/])+)\*\/')      
            param = re.compile(r'(?:if\s*\(\!)?\s*sf_get'
                               '(?P<type>bool|int|float)'
                               '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                               '\s*\&(?P<var>[\w\_\[\]]+)\s*[\)]\s*[\)]?\s*'
                               '(?:[\{]|' # either \{ or
                               '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|' 
                               'sf_[^\;]+)?' # or sf_error
                               '[\;])\s*' # ending with ;
                               '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                               '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
            params = re.compile(r'sf_get(?P<type>bools|ints|floats|strings)'
                                '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                                '\s*(?P<var>[\w\_\[\]]+)\s*\,'
                                '\s*(?P<size>[\w\_]+)\s*\)\s*'
                                '[^\;\{]*[\;\{]\s*' # ending with ; or {
                                '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                                '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
            param2 = re.compile(r'sf_get(?P<type>bool|int|float|string)\s*'
                                '\([^/]+\/\*\<\s*(?P<name>[\w\#]+)'
                                '(?:=(?P<default>\S+))?'
                                '\s*(?P<desc>[^\>]+)\>\*\/')
            params2 = re.compile(r'sf_get(?P<type>bools|ints|floats|strings)'
                                 '\s*\([^\,]+\,[^\,]+\,'
                                 '\s*(?P<size>[\w\_]+)\s*\)[^/]+'             
                                 '\/\*\<\s*(?P<name>[\w\#]+)'
                                 '(?:=(?P<default>\S+))?'
                                 '\s*(?P<desc>[^\>]+)\>\*\/')
            stringpar = re.compile(r'sf_getstring\s*\(\s*\"(?P<name>\w+)\"'
                                   '[^\;\{]*[\;\{]\s*(?:\/\*'
                                   '\s*(?P<desc>(?:[^*]|\*[^/])+)\*\/)?')
            synopsis = re.compile(r'\s*Takes\s*\:\s*((?:[^\n]|[\n][^\n])+)'
                                  '((?:.|\n)*)$')
            inpout = re.compile(r'\s*(?P<name>\w+)\s*=\s*'
                                'sf_(?P<io>input|output)'
                                '\s*\(\s*\"(?P<tag>\w+)\"')            
            version = re.compile(r'\/\*\s*\$Id\:\s*(.+\S)\s*\$\s*\*\/')
    name = rsfprefix + re.sub('^M','',os.path.basename(file))
    if lang[0] =='f':
        name = re.sub('\.f\d*$','',name)
    else:
        name = re.sub('\.cc?$','',name)
    src = open(file,"r")   # open source
    text = string.join(src.readlines(),'')
    src.close()
    first = comment.match(text)
    if first:
        tops = string.split(first.group(1),"\n")
        desc = string.lstrip(tops.pop(0))
        first = string.join(tops,"\n")
    else:
        desc = None
    prog = rsfprog(name,file,desc)
    file = re.sub('^[^\/]*\/','',file)
    out.write("%s = rsfdoc.rsfprog('%s','%s','''%s''')\n" %
              (name,name,file,desc))
    parline = ''
    if params:
        pars = params.findall(text)
    else:
        pars = []
    for par in pars:
        type = par[0]
        parname = par[1]
        size = par[3]
        range = par[4]
        desc = par[5] + ' [%s]' % size

        prog.par(parname,rsfpar(type,None,range,desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                  (name,parname,type,'',range,desc))
        parline = parline + " %s=" % parname
    pars = param.findall(text)
    for par in pars:
        if lang == 'f90':
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

        prog.par(parname,rsfpar(type,default,range,desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                  (name,parname,type,default,range,desc))
        parline = parline + " %s=%s" % (parname,default)
    if params2:
        pars = params2.findall(text)
    else:
        pars = []
    for par in pars:
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
    if param2: 
        pars = param2.findall(text)
    else:
        pars = []
    for par in pars:
        type = par[0]
        parname = par[1]
        default = par[2]
        range = '' # for now
        desc = par[3]

        prog.par(parname,rsfpar(type,default,range,desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                  (name,parname,type,default,range,desc))
        parline = parline + " %s=%s" % (parname,default)
    if stringpar:
        pars = stringpar.findall(text)
    else:
        pars = []
    for par in pars:
        type = 'string '
        parname = par[0]
        desc = par[1]
        prog.par(parname,rsfpar(type,desc=desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('%s',desc='''%s'''))\n" %
                  (name,parname,type,desc))
        parline = parline + " %s=" % (parname)
    files = inpout.findall(text)
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
            snps = snps + iostring
        valid[filename]=0
    if re.match(rsfplotprefix,name):
        snps = snps + ' > plot.' + rsfplotsuffix
    snps = snps + parline
    vers = version.search(text)
    if vers:
        prog.version(vers.group(1))
        out.write("%s.version('''%s''')\n" % (name,vers.group(1)))
    if not first:
        first = ''
    else:
        info = synopsis.match(first)
        if info:
            snps = snps + ' ' + string.lstrip(info.group(1))
            first = string.lstrip(info.group(2))
    prog.synopsis(snps,first)
    out.write("%s.synopsis('''%s''','''%s''')\n" % (name,snps,first))
    out.write("rsfdoc.progs['%s']=%s\n\n" % (name,name))

def cli(rsfprefix = 'sf',rsfplotprefix='vp'):
    import getopt
    import rsfprog

    this = sys.argv.pop(0)
    class BadUsage: pass

    try:
        opts, args = getopt.getopt(sys.argv, 'k:w:m:l:r:')
        hdir = None
        mdir = None
        ldir = None
        rep = os.getcwd()
        for opt, val in opts:
            if opt == '-w':
                hdir = val
            if opt == '-m':
                mdir = val
            if opt == '-l':
                ldir = val
            if opt == '-r':
                rep = val
            if opt == '-k':
                val = val.lower()
                doc = ''
                for prog in progs.keys():
                    desc = progs[prog].desc
                    if re.search(val,desc.lower()):
                        doc = doc + "%s: %s\n" % (bold(prog),desc)
                pydoc.pager(doc)
                return
    
        if not args:
            if hdir:
                html(hdir)
                for prog in progs.keys():
                    main = progs.get(prog)
                    if main:
                        main.html(hdir,rep)
            else:
                raise BadUsage

        for prog in args:
            if     not re.match(rsfprefix,prog) and \
                   not re.match(rsfplotprefix,prog):
                prog = rsfprefix + prog
            main = progs.get(prog)
            if main:
                if hdir:
                    main.html(hdir,rep)
                elif mdir:
                    main.wiki(mdir)
                elif ldir:
                    main.latex(ldir)
                else:
                    main.document()
            else:
                print "No program %s in RSF." % prog

    except (getopt.error, BadUsage):
        print '''sfdoc - the RSF documentation tool
        
%(prog)s <prog1> <prog2> ... 
    Show documentation on programs.

%(prog)s -w <dir> [-r <rep>] <prog1> <prog2> ... 
    Write program HTML documentaton in <dir> directory, optional <rep> referes to repository.

%(prog)s -m <dir>  <prog1> <prog2> ...
    Write program documentaton in MediaWiki format in <dir> directory.

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

# 	$Id$
