import os, stat, sys, types, commands, re, string, urllib
import rsfdoc
import rsfprog
import rsfconf

##############################################################################
# BEGIN STANDARD SCons SCRIPT HEADER
#
# This is the cut-and-paste logic so that a self-contained script can
# interoperate correctly with different SCons versions and installation
# locations for the engine.  If you modify anything in this section, you
# should also change other scripts that use this same header.
##############################################################################

# Strip the script directory from sys.path() so on case-insensitive
# (WIN32) systems Python doesn't think that the "scons" script is the
# "SCons" package.  Replace it with our own library directories
# (version-specific first, in case they installed by hand there,
# followed by generic) so we pick up the right version of the build
# engine modules if they're in either directory.

script_dir = sys.path[0]

if script_dir in sys.path:
    sys.path.remove(script_dir)

libs = []

if os.environ.has_key("SCONS_LIB_DIR"):
    libs.append(os.environ["SCONS_LIB_DIR"])

prefs = []

if sys.platform == 'win32':
    # sys.prefix is (likely) C:\Python*;
    # check only C:\Python*.
    prefs.append(sys.prefix)
else:
    # On other (POSIX) platforms, things are more complicated due to
    # the variety of path names and library locations.  Try to be smart
    # about it.
    if script_dir == 'bin':
        # script_dir is `pwd`/bin;
        # check `pwd`/lib/scons*.
        prefs.append(os.getcwd())
    else:
        if script_dir == '.' or script_dir == '':
            script_dir = os.getcwd()
        head, tail = os.path.split(script_dir)
        if tail == "bin":
            # script_dir is /foo/bin;
            # check /foo/lib/scons*.
            prefs.append(head)

    head, tail = os.path.split(sys.prefix)
    if tail == "usr":
        # sys.prefix is /foo/usr;
        # check /foo/usr/lib/scons* first,
        # then /foo/usr/local/lib/scons*.
        prefs.append(sys.prefix)
        prefs.append(os.path.join(sys.prefix, "local"))
    elif tail == "local":
        h, t = os.path.split(head)
        if t == "usr":
            # sys.prefix is /foo/usr/local;
            # check /foo/usr/local/lib/scons* first,
            # then /foo/usr/lib/scons*.
            prefs.append(sys.prefix)
            prefs.append(head)
        else:
            # sys.prefix is /foo/local;
            # check only /foo/local/lib/scons*.
            prefs.append(sys.prefix)
    else:
        # sys.prefix is /foo (ends in neither /usr or /local);
        # check only /foo/lib/scons*.
        prefs.append(sys.prefix)

    prefs = map(lambda x: os.path.join(x, 'lib'), prefs)

libs.extend(map(lambda x: os.path.join(x, 'scons'), prefs))

sys.path = libs + sys.path

##############################################################################
# END STANDARD SCons SCRIPT HEADER
##############################################################################

from SCons.Environment import Environment
from SCons.Util import WhereIs
from SCons.Builder import Builder
from SCons.Action import Action
from SCons.Scanner import Base
from SCons.Options import Options
from SCons.Node.FS import default_fs 

##############################################################################
# BEGIN CONFIGURATION VARIABLES
##############################################################################

#prefix for sf commands
sfprefix = 'sf'
#suffix for rsf files
sfsuffix = '.rsf'
# suffix for vplot files
vpsuffix = '.vpl'
# suffix for eps files
pssuffix = '.eps'
# path bor binary files
datapath = os.environ.get('DATAPATH')
if not datapath:
    try:
        file = open('.datapath','r')
    except:
        try:
            file = open(os.path.join(os.environ.get('HOME'),'.datapath'),'r')
        except:
            file = None
    if file:
        for line in file.readlines():
            check = re.match("(?:%s\s+)?datapath=(\S+)" % os.uname()[1],line)
            if check:
                datapath = check.group(1)
        file.close()
    if not datapath:
        datapath = os.path.join(os.environ.get('HOME'),'')

# directory tree for executable files
top = os.environ.get('RSFROOT')
bindir = os.path.join(top,'bin')
libdir = os.path.join(top,'lib')
incdir = os.path.join(top,'include')

resdir = './Fig'

latex = None
bibtex = None
rerun = None

# temporary (I hope)
sep = os.path.join(os.environ.get('SEP'),'bin/')

##############################################################################
# END CONFIGURATION VARIABLES
##############################################################################

def collect_exe(dir):
    "Make a list of executable files in a directory"
    def isexe(file):
        "Check if a file is executable" # Posix only?
        return (os.stat(file)[stat.ST_MODE] &
                (stat.S_IXUSR|stat.S_IXGRP|stat.S_IXOTH))
    exe = []
    for file in os.listdir(dir):
        file = os.path.join(dir,file)
        if os.path.isfile(file) and isexe(file):
            exe.append(file)
    return exe

#############################################################################
# CUSTOM BUILDERS
#############################################################################

#def clean(target=None,source=None,env=None):
#    for junk in env['junk']:
#        if (os.path.isfile (junk)):
#            try:
#                os.unlink(junk)
#            except:
#                pass
#    return 0

def silent(target=None,source=None,env=None):
    return None

ppi = 72 # points per inch resolution
def pstexpen(target=None,source=None,env=None):
    "Convert vplot to EPS"
    vplot = str(source[0])
    eps = str(target[0])
    space = os.environ.get('PSBORDER')
    if not space:
        space=0.
    else:
        space=float(space)
    opts = os.environ.get('PSTEXPENOPTS')
    if not opts:
        opts = ''
    pstexpenopts = env.get('opts')
    if not pstexpenopts:
        pstexpenopts = 'color=n fat=1 fatmult=1.5 invras=y'
    opts = ' '.join([opts,pstexpenopts])
    print opts
    head = string.split(
        commands.getoutput(sep +
                           "vppen big=n stat=l %s < %s | %s -1" %
                           (opts,vplot,WhereIs('head'))))
    bb = []
    for x in (7, 12, 9, 14):
        bb.append(int((float(head[x])-space)*ppi))
    try:
        file = open(eps,"w")
        file.write("%\!PS-Adobe-2.0 EPSF-2.0\n")
        file.write("%%%%BoundingBox: %d %d %d %d\n" % tuple(bb))
        file.write(commands.getoutput(sep + "pspen size=a tex=y %s < %s" %
                                      (opts,vplot)))
        file.write("\n")
        file.close()
    except:
        return 1
    return 0

# should not LATEX be in env?
#if WhereIs('dvips'):
#    latex = WhereIs('latex')
#else:
# latex = WhereIs('pdflatex')
bibtex = WhereIs('bibtex')
rerun = re.compile(r'\bRerun')

def latex2dvi(target=None,source=None,env=None):
    "Convert LaTeX to DVI/PDF"
    latex = env.get('latex',WhereIs('pdflatex'))
    tex = str(source[0])
    dvi = str(target[0])
    stem = re.sub('\.[^\.]+$','',dvi)    
    run = string.join([latex,tex],' ')
    # First latex run
    if os.system(run):
        return 1
    # Check if bibtex is needed
    aux = open(stem + '.aux',"r")    
    for line in aux.readlines():
        if re.search("bibdata",line):
            os.system(string.join([bibtex,stem],' '))
            os.system(run)
            os.system(run)
            break        
    aux.close()
    # (Add makeindex later)
    # Check if rerun is needed
    for i in range(3): # repeat 3 times at most
        done = 1
        log = open(stem + '.log',"r")
        for line in log.readlines():
            if rerun.search(line):
                done = 0
                break
        log.close()
        if done:
            break
        os.system(run)
    return 0

def retrieve(target=None,source=None,env=None):
    "Fetch data from the web"
    dir = env['dir']
    for file in map(str,target):
        urllib.urlretrieve(string.join([dir,file],os.sep),file)
    return 0

View = Builder(action = sep + "xtpen $SOURCES",src_suffix=vpsuffix)
# Klean = Builder(action = Action(clean,silent,['junk']))
Build = Builder(action = Action(pstexpen,varlist=['opts']),
                src_suffix=vpsuffix,suffix=pssuffix)
epstopdf = WhereIs('epstopdf')
if epstopdf:
    PDFBuild = Builder(action = epstopdf + " $SOURCES",
		       src_suffix=pssuffix,suffix='.pdf')
Retrieve = Builder(action = Action(retrieve,varlist=['dir']))

#if WhereIs('dvips'):
#    dvips = 1
#    Dvi = Builder(action = Action(latex2dvi),
#              src_suffix=['.tex','.ltx'],suffix='.dvi')
#    Ps = Builder(action = WhereIs('dvips') + " -Ppdf -G0 -o $TARGET $SOURCES",
#                 src_suffix='.dvi',suffix='.ps')
#    Pdf = Builder(action = WhereIs('ps2pdf') + " $SOURCES",
#                  src_suffix='.ps',suffix='.pdf')
#    ressuffix = '.ps'
#else:
dvips = 0
Pdf = Builder(action = Action(latex2dvi,varlist=['latex']),
              src_suffix=['.tex','.ltx'],suffix='.pdf')

acroread = WhereIs('acroread')
if acroread:
    Read = Builder(action = acroread + " $SOURCES",src_suffix='.pdf')
ressuffix = '.pdf'

#############################################################################
# CUSTOM SCANNERS
#############################################################################

isplot = None
def getplots(node,env,path):
    global isplot, ressufix
    if not isplot:
        isplot = re.compile(r'\\(?:side)?plot\s*\{([^\}]+)')
    contents = node.get_contents()
    plots = isplot.findall(contents)
    return map(lambda x: os.path.join(resdir,x) + ressuffix,plots)

Plots = Base(name='Plots',function=getplots,skeys=['.tex','.ltx'])

#############################################################################
# PLOTTING COMMANDS
#############################################################################

combine ={
    'SideBySideAniso': lambda n:
    sep + "vppen yscale=%d vpstyle=n gridnum=%d,1 $SOURCES" % (n,n),
    'OverUnderAniso': lambda n:
    sep + "vppen xscale=%d vpstyle=n gridnum=1,%d $SOURCES" % (n,n),
    'SideBySideIso': lambda n:
    sep + "vppen size=r vpstyle=n gridnum=%d,1 $SOURCES" % n,
    'OverUnderIso': lambda n:
    sep + "vppen size=r vpstyle=n gridnum=1,%d $SOURCES" % n,
    'TwoByTwo': lambda n:
    sep + "vppen size=r vpstyle=n gridnum=2,2 $SOURCES",
    'Overlay': lambda n:
    sep + "vppen erase=o vpstyle=n $SOURCES",
    'Movie': lambda n:
    sep + "vppen vpstyle=n $SOURCES"
    }

#############################################################################

display = os.environ.get('DISPLAY')
if display:
    host = re.sub(':[\d\.]*$','',display)
    if host == '':
        host = 'localhost'
    os.system('xhost ' + host)

class Project(Environment):
    def __init__(self,**kw):
        apply(Environment.__init__,(self,),kw)
        # Add f90 later
        opts = Options(os.path.join(libdir,'rsfconfig.py'))
        rsfconf.options(opts)
        opts.Update(self)
        dir = os.path.basename(os.getcwd())
        self.path = datapath + dir + os.sep
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        self.Append(ENV={'DATAPATH':self.path,
                         'DISPLAY':display,
                         'RSFROOT':top},
                    BUILDERS={'View':View,
#                              'Clean':Klean,
                              'Build':Build,
#                              'PDFBuild':PDFBuild,
                              #                              'Dvi':Dvi,
                              #                              'Ps':Ps,
#                              'Read':Read,
                              'Retrieve':Retrieve},
                    SCANNERS=[Plots],
                    LIBPATH=[libdir],
                    CPPPATH=[incdir],
                    LIBS=['rsf','m'],
                    PROGSUFFIX='.x')
	if acroread:
	    self.Append(BUILDERS={'Read':Read})
	if epstopdf:
	    self.Append(BUILDERS={'PDFBuild':PDFBuild})
        self['PROGPREFIX']=''
        self.view = []
        self.figs = []
        self.pdfs = []
    def Flow(self,target,source,flow,stdout=1,stdin=1,
             suffix=sfsuffix,prefix=sfprefix,src_suffix=sfsuffix):
        if not flow:
            return None        
        sources = []
        if source:
            if type(source) is types.ListType:
                files = source
            else:
                files = string.split(source)
            for file in files:
                if ('.' not in file):
                    file = file + src_suffix
                sources.append(file)
        else:
            stdin=0
        lines = string.split(flow,';')
        steps = []
        for line in lines:
            substeps = []
            sublines = string.split(line,'|')
            for subline in sublines:           
                pars = string.split(subline)
                # command is assumed to be always first in line
                command = pars.pop(0)
                # check if this command is in our list
                rsfprog = prefix + command            
                if rsfdoc.progs.has_key(rsfprog):
                    command = os.path.join(bindir,rsfprog)
                    sources.append(command)
                elif re.match(r'[^/]+\.x$',command): # local program
                    command = os.path.join('.',command)
                #<- check for par files and add to the sources
                for par in pars:
                    if re.match("^par=",par):
                        sources.append(default_fs.File(par[4:]))
                #<<- assemble the command line
                pars.insert(0,command)
                substeps.append(string.join(pars,' '))
            #<-
            steps.append(string.join(substeps," | "))
        #<- assemble the pipeline
        command = string.join(steps," ;\n")
        if stdout==1:
            command = command + " > $TARGET"
        elif stdout==0:
            command = command + " >/dev/null"
        if stdin:
            command = "< $SOURCE " + command
        targets = []
        if type(target) is types.ListType:
            files = target
        else:
            files = string.split(target)
        for file in files:
            if not re.search(suffix + '$',file):
                file = file + suffix
            targets.append(file)
        if suffix == sfsuffix:            
            datafiles = [] 
            for target in targets:
                if os.sep not in target:
                    datafile = self.path + target + '@'
                    datafiles.append(datafile)
            targets = targets + datafiles
        return self.Command(targets,sources,command)
    def Plot (self,target,source,flow,suffix=vpsuffix,**kw):
        kw.update({'suffix':suffix})
        return apply(self.Flow,(target,source,flow),kw)
    def Result(self,target,source,flow,suffix=vpsuffix,
               pstexpen=None,**kw):
        target2 = os.path.join(resdir,target)
        if flow:
#            kw.update({'clean':clean,'suffix':suffix})
            plot = apply(self.Plot,(target2,source,flow),kw)
            self.Default (plot)
            self.view.append(self.View(target + '.view',plot))
            build = self.Build(target2 + pssuffix,plot,opts=pstexpen)
            self.figs.append(build)
            self.Alias(target + '.build',build)
        else:
            plot = None
            build = target2 + pssuffix
	if epstopdf:
	    buildPDF = self.PDFBuild(target2,build)
	    self.pdfs.append(buildPDF)
	    self.Alias(target + '.buildPDF',buildPDF)
        return plot
    def Combine(self,target,source,how,result=0,vppen=None,**kw):
        if not type(source) is types.ListType:
            source = string.split(source)
        flow = apply(combine[how],[len(source)])
        if vppen:
            flow = flow + ' ' + vppen
        kw.update({'src_suffix':vpsuffix,'stdin':0})
        if result:
            return apply(self.Result,(target,source,flow),kw)
        else:
            return apply(self.Plot,(target,source,flow),kw)
    def End(self):
        self.Alias('view',self.view)
#        self.Alias('clean',self.Clean('clean',None,junk=self.junk))
        if self.figs: # if any results
            build = self.Alias('build',self.figs)
        if self.pdfs:
            buildPDF = self.Alias('buildPDF',self.pdfs)
        self.Append(BUILDERS={'Pdf':Pdf})
        if os.path.isfile('paper.tex'): # if there is a paper
            if dvips:
                self.paper = self.Dvi(target='paper',source='paper.tex')
                self.Alias('dvi',self.paper)
                self.Alias('ps',self.Ps('paper'))
                self.Alias('pdf',self.Pdf('paper'))
            else:
                self.paper = self.Pdf(target='paper',source='paper.tex')
                self.Alias('pdf',self.paper)
            self.paper.target_scanner = Plots
	    if acroread:
		self.Alias('read',self.Read('paper'))
    def Fetch(self,file,dir):
        return self.Retrieve(file,None,dir=dir)

# Default project
project = Project()
def Flow(target,source,flow,**kw):
    return apply(project.Flow,(target,source,flow),kw)
def Plot (target,source,flow,**kw):
    return apply(project.Plot,(target,source,flow),kw)
def Result(target,source,flow,**kw):
    return apply(project.Result,(target,source,flow),kw)
def Combine(target,source,how,**kw):
    return apply(project.Combine,(target,source,how),kw)
def Fetch(file,dir):
    return project.Fetch(file,dir)
def End():
    project.End()

if __name__ == "__main__":
     import pydoc
     pydoc.help(Project)
     
# 	$Id: rsfproj.py,v 1.22 2004/03/27 03:28:58 fomels Exp $	
