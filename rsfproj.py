import os, stat, sys, types, commands, re, string, urllib
import rsfdoc
import rsfprog
import rsfconfig

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
from SCons.Script.SConscript import Default, Clean
from SCons.Util import WhereIs
from SCons.Builder import Builder
from SCons.Action import Action
from SCons.Scanner import Base

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
# path bor binary files (change later for compliance with SEPlib)
datapath = os.path.join(os.path.join(os.environ.get('HOME'),'scr'),'')

# directory tree for executable files
top = os.environ.get('RSFROOT')
bindir = os.path.join(top,'bin')

resdir = './Fig'

latex = None
bibtex = None
rerun = None

# temporary (I hope)
sep = os.path.join(top,"../SEP/bin/")

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

def clean(target=None,source=None,env=None):
    for junk in env['junk']:
        if (os.path.isfile (junk)):
            try:
                os.unlink(junk)
            except:
                pass
    return 0

def silent(target=None,source=None,env=None):
    return None

ppi = 72 # points per inch resolution
def pstexpen(target=None,source=None,env=None):
    "Convert vplot to EPS"
    vplot = str(source[0])
    eps = str(target[0])
    space = os.environ.get('PSBORDER')
    if (not space):
        space=0.
    else:
        space=float(space)
    opts = os.environ.get('PSTEXPENOPTS')
    if (not opts):
        opts = ''
    opts = opts + "color=n fat=1 fatmult=1.5 invras=y"
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
latex = WhereIs('pdflatex')
bibtex = WhereIs('bibtex')
rerun = re.compile(r'\bRerun')

def latex2dvi(target=None,source=None,env=None):
    "Convert LaTeX to DVI"
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
        urllib.urlretrieve(string.join([dir,file],'/'),file)
    return 0

View = Builder(action = sep + "xtpen $SOURCES",src_suffix=vpsuffix)
Klean = Builder(action = Action(clean,silent,['junk']))
Build = Builder(action = Action(pstexpen),src_suffix=vpsuffix,suffix=pssuffix)
PDFBuild = Builder(action = WhereIs('epstopdf') + " $SOURCES",
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
Pdf = Builder(action = Action(latex2dvi),
              src_suffix=['.tex','.ltx'],suffix='.pdf')

Read = Builder(action = WhereIs('acroread') + " $SOURCES",
               src_suffix='.pdf')
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
    'Overlay': lambda n:
    sep + "vppen erase=o vpstyle=n $SOURCES"
    }

#############################################################################

class Project(Environment):
    def __init__(self,**kw):
        global view, dvips
        Environment.__init__(self,**kw)
        self.Append(ENV={'DATAPATH':datapath,
                         'DISPLAY':os.environ.get('DISPLAY')},
                    BUILDERS={'View':View,
                              'Clean':Klean,
                              'Build':Build,
                              'PDFBuild':PDFBuild,
#                              'Dvi':Dvi,
#                              'Ps':Ps,
                              'Pdf':Pdf,
                              'Read':Read,
                              'Retrieve':Retrieve},
                    SCANNERS=[Plots])
# Add f90 later
        self.view = []
        self.figs = []
        self.pdfs = []
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
            self.Alias('read',self.Read('paper'))
            self.junk = ['paper.aux','paper.log','paper.bbl',
                         'paper.blg','paper.ps','paper.dvi']
            Clean('paper.pdf',self.junk)
        else:
            self.junk = []
    def Flow(self,target,source,flow,clean=1,stdout=1,stdin=1,
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
                rsfprog = prefix + command            
                # check if this command is in our list
                if rsfdoc.progs.has_key(rsfprog):
                    command = os.path.join(bindir,rsfprog)
                    sources.append(command)
                #<- check for par files and add to the sources
                for par in pars:
                    if re.match("^par=",par):
                        sources.append(par[4:])
                #<<- assemble the command line
                pars.insert(0,command)
                substeps.append(string.join(pars,' '))
            #<-
            steps.append(string.join(substeps," | "))
        #<- assemble the pipeline
        command = string.join(steps," ;\n")
        if stdout:
            command = command + " > $TARGET"
        else:
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
        if clean:
            self.junk.append(target)
        return self.Command(targets,sources,command)
    def Plot (self,target,source,flow,suffix=vpsuffix,**kw):
        return self.Flow(target,source,flow,suffix=suffix,**kw)
    def Result(self,target,source,flow,clean=0,suffix=vpsuffix,**kw):
        target2 = os.path.join(resdir,target)
        if flow:
            plot = self.Plot(target2,source,flow,
                             clean=clean,suffix=suffix,**kw)
            Default (plot)
            self.view.append(self.View(target + '.view',plot))
            build = self.Build(target2 + pssuffix,plot)
            self.figs.append(build)
            self.Alias(target + '.build',build)
        else:
            plot = None
            build = target2 + pssuffix
        buildPDF = self.PDFBuild(target2,build)
        self.pdfs.append(buildPDF)
        self.Alias(target + '.buildPDF',buildPDF)
        return plot
    def Combine(self,target,source,how,result=0):
        if not type(source) is types.ListType:
            source = string.split(source)
        flow = apply(combine[how],[len(source)])
        if result:
            return self.Result(target,source,flow,src_suffix=vpsuffix,stdin=0)
        else:
            return self.Plot(target,source,flow,src_suffix=vpsuffix,stdin=0)  
    def End(self):
        self.Alias('view',self.view)
        self.Alias('clean',self.Clean('clean',None,junk=self.junk))
        if self.figs: # if any results
            build = self.Alias('build',self.figs)
        if self.pdfs:
            buildPDF = self.Alias('buildPDF',self.pdfs)
    def Fetch(self,file,dir):
        return self.Retrieve(file,None,dir=dir)

if __name__ == "__main__":
     proj = Project()
     import pydoc
     pydoc.help(Project)
     
