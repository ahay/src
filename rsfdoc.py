import pydoc
import re, sys, os, string

progs = {}

def selfdoc(target=None,source=None,env=None):
    src = str(source[0])
    doc = open(str(target[0]),"w")
    getprog(src,doc)
    doc.close()

def bold(text):
    """Format a string in bold by overstriking."""
    return string.join(map(lambda ch: ch + "\b" + ch, text),'')

def underline(text):
    """Format a string in underline by overstriking."""
    return string.join(map(lambda ch: ch + "\b_", text),'')

def section(head,body):
    text = string.join(map(lambda line: "\t" + line,
                           string.split(body,"\n")),"\n")
    return bold(string.upper(head)) + "\n" + text + "\n"

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
        self.pars = {}
    def synopsis (self,snps,cmts):
        self.snps = snps
        self.cmts = cmts
    def par (self,name,par):
        self.pars[name] = par
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
        doc = doc + section('source',self.file)
        pydoc.pager(doc)

comment = None
param = None
stringpar = None
synopsis = None

rsfprefix = 'sf'

def getprog(file,out):
    global comment, param, synopsis, stringpar
    if not comment:
        comment = re.compile(r'\/\*((?:[^*]|\*[^/])+)\*\/')
        param = re.compile(r'(?:if\s*\(\!)?sf_get(?P<type>bool|int|float)\s*'
                           '\(\s*\"(?P<name>\w+)\"\s*\,'
                           '\s*\&(?P<var>[\w\[\]]+)\s*'
                           '\)\s*[\)]?\s*[\{]?\s*'
                           '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|'
                           'sf_error[^\;]+)?'
                           '\;\s*(?:\/\*\s*(?P<desc>(?:[^*]|\*[^/])+)\*\/)?')
        stringpar = re.compile(r'sf_getstring\s*\(\s*\"(?P<name>\w+)\"[^\;]*\;'
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
    out.write("%s = rsfdoc.rsfprog('%s','%s','%s')\n" %
              (name,name,file,desc))
    pars = param.findall(text)
    parline = ''
    for par in pars:
        type = par[0]
        parname = par[1]
        default = par[3]
        desc = par[4]
        range = ''
        if (type == 'bool'):
            if (default == 'true'):
                default = 'y'
            else:
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

if __name__ == "__main__":
    junk = open('junk.py',"w")
    junk.write("import rsfdoc\n\n")
    junk.write("rsfprog = {}\n")
    getprog('seis/main/dd.c',junk)
    junk.write("sfdd.document()\n\n")
    junk.write("print sfdd.prog\n\n")
    junk.close()
    #
    import junk
    #
    os.unlink("junk.py")
    os.unlink("junk.pyc")


