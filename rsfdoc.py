import pydoc
import re, sys, os, string, glob, commands

progs = {}

def subdirs():
    return filter(os.path.isdir,glob.glob('[a-z]*'))

def use(target=None,source=None,env=None):
    doc = open(str(target[0]),'w')
    doc.write('import rsfdoc\n\n')
    os.chdir('book')
    for book in subdirs():
        os.chdir(book)
        for chapter in subdirs():
            os.chdir(chapter)
            for project in subdirs():
                os.chdir(project)
                (status,progs) = commands.getstatusoutput('scons -s .sf_uses')
                if status:
                    print ('No uses found in book/%s/%s/%s/: %s' %
                           (book,chapter,project,progs))
                else:
                    for prog in progs.split():
                        doc.write('rsfdoc.progs["%s"].use("%s","%s","%s")\n' %
                                  (prog,book,chapter,project))
                os.chdir('..')
            os.chdir('..')
        os.chdir('..')
    os.chdir('..')
    doc.close()

def selfdoc(target=None,source=None,env=None):
    src = str(source[0])
    doc = open(str(target[0]),"w")
    rsfprefix = env.get('rsfprefix','sf')
    getprog(src,doc,rsfprefix)
    doc.close()

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

class rsfpar:
    def __init__(self,type,default='',range='',desc=''):
        self.type = underline(type) + " "
        self.default = "=" + str(default)
        self.range = " " + range
        self.desc = "\t" + desc + "\n"
    def show(self,name):
        return self.type + bold(name + self.default) + self.range + self.desc

class rsfprog:
    def __init__(self,name,file,desc=None):
        self.name = name
        self.file = file
        self.desc = desc
        self.snps = None
        self.cmts = None
        self.also = None
        self.uses = {}
        self.pars = {}
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
        pydoc.pager(doc)
    def html(self,dir):
        file = open (os.path.join(dir,self.name + '.html'),'w')
        file.write(page(self.name,'Test'))
        file.close()

comment = None
param = None
stringpar = None
synopsis = None

def getprog(file,out,rsfprefix = 'sf'):
    global comment, param, synopsis, stringpar
    if not comment:
        comment = re.compile(r'\/\*((?:[^*]|\*[^/])+)\*\/')
        param = re.compile(r'(?:if\s*\(\!)?sf_get(?P<type>bool|int|float)\s*'
                           '\(\s*\"(?P<name>\w+)\"\s*\,'
                           '\s*\&(?P<var>[\w\_\[\]]+)\s*[\)]\s*[\)]?\s*'
                           '(?:[\{]|' # either \{ or
                           '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|' # var=def
                           'sf_[^\;]+)?' # or sf_error
                           '[\;])\s*' # ending with ;
                           '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                           '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
        stringpar = re.compile(r'sf_getstring\s*\(\s*\"(?P<name>\w+)\"'
                               '[^\;\{]*[\;\{]'
                               '\s*(?:\/\*\s*(?P<desc>(?:[^*]|\*[^/])+)\*\/)?')
        synopsis = re.compile(r'\s*Takes\s*\:\s*((?:[^\n]|[\n][^\n])+)'
                              '((?:.|\n)*)$')
    name = rsfprefix + re.sub('^M','',os.path.basename(file))
    name = re.sub('.c$','',name)
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
    out.write("%s = rsfdoc.rsfprog('%s','%s','%s')\n" %
              (name,name,file,desc))
    pars = param.findall(text)
    parline = ''
    for par in pars:
        type = par[0]
        parname = par[1]
        default = par[3]
        range = par[4]
        desc = par[5]
        if (type == 'bool'):
            if (default == 'true'):
                default = 'y'
            elif (default == 'false'):
                default = 'n'
            type = 'bool  ' # to align with string
            range = '[y/n]'
        elif (type == 'int'):
            type = 'int   ' # to align with string
        elif (type == 'float'):
            type = 'float ' # to align with string
        prog.par(parname,rsfpar(type,default,range,desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('%s','%s','%s','''%s'''))\n" %
                  (name,parname,type,default,range,desc))
        parline = parline + " %s=%s" % (parname,default)
    pars = stringpar.findall(text)
    for par in pars:
        type = 'string'
        parname = par[0]
        desc = par[1]
        prog.par(parname,rsfpar("string",desc=desc))
        out.write("%s.par('%s',rsfdoc.rsfpar('string',desc='''%s'''))\n" %
                  (name,parname,desc))
        parline = parline + " %s=" % (parname)
    if first:
        info = synopsis.match(first)
        if info:
            snps = name + " " + string.lstrip(info.group(1)) + parline
            cmts = string.lstrip(info.group(2))
            prog.synopsis(snps,cmts)
            out.write("%s.synopsis('''%s''','''%s''')\n" % (name,snps,cmts))
    out.write("rsfdoc.progs['%s']=%s\n\n" % (name,name))


def cli(rsfprefix = 'sf'):
    import getopt
    import rsfprog

    this = sys.argv.pop(0)
    class BadUsage: pass

    try:
        opts, args = getopt.getopt(sys.argv, 'k:w:')
        dir = None
        for opt, val in opts:
            if opt == '-w':
                dir = val
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
            raise BadUsage

        for prog in args:
            if not re.match(rsfprefix,prog):
                prog = rsfprefix + prog
            main = progs.get(prog)
            if main:
                if dir:
                    main.html(dir)
                else:
                    main.document()
            else:
                print "No program %s in RSF." % prog

    except (getopt.error, BadUsage):
        print '''sfdoc - the RSF documentation tool
        
%(prog)s <prog1> <prog2> ... 
    Show documentation on programs.

%(prog)s <prog1> <prog2> ... -w <dir>
    Write program documentaton in <dir> directory.

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

# 	$Id: rsfdoc.py,v 1.11 2004/03/31 03:16:33 fomels Exp $	
