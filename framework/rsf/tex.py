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

from __future__ import print_function, division, absolute_import
import os, re, glob, string, types, pwd, shutil
import io, token, tokenize, cgi, sys, keyword

import rsf.conf, rsf.path, rsf.latex2wiki
import rsf.doc
import rsf.prog

import SCons

# The following adds all SCons SConscript API to the globals of this module.
version = list(map(int,SCons.__version__.split('.')[:3]))
if version[0] >= 1 or version[1] >= 97 or (version[1] == 96 and version[2] >= 90):
    from SCons.Script import *
else:
    import SCons.Script.SConscript
    globals().update(SCons.Script.SConscript.BuildDefaultGlobals())

SCons.Defaults.DefaultEnvironment(tools = [])

#############################################################################
# CONFIGURATION VARIABLES
#############################################################################

# To do: the configuration process should know about these dependencies!

bibtex      = WhereIs('bibtex')
makeindex   = WhereIs('makeindex')
acroread    = WhereIs('acroread')
pdfread     = acroread or WhereIs('kpdf') or WhereIs('evince') or \
    WhereIs('xpdf') or WhereIs('gv') or WhereIs('open')
pdftops     = WhereIs('pdftops')
epstopdf    = WhereIs('epstopdf') or WhereIs('a2ping')
if epstopdf:
    latex       = WhereIs('pdflatex')
    ressuffix = '.pdf'
else:
    latex       = WhereIs('latex')
    ressuffix = '.eps'
fig2dev     = WhereIs('fig2dev')
latex2html  = WhereIs('latex2html')
pdf2ps      = WhereIs('pdf2ps')
ps2eps      = WhereIs('gs') or WhereIs('ps2epsi')
pstoimg     = WhereIs('pstoimg')
mathematica = WhereIs('mathematica')
if mathematica:
    mathematica = WhereIs('math')
matlab      = WhereIs('matlab')
gnuplot     = WhereIs('gnuplot') or WhereIs('gnuplot44')
sage        = WhereIs('sage')

try:
    import numpy, pylab
    haspylab=1
except:
    haspylab=0

vpsuffix  = '.vpl'
pssuffix  = '.eps'
itype = os.environ.get('IMAGE_TYPE','png')

rerun = re.compile(r'\bRerun')

#############################################################################
# REGULAR EXPRESSIONS
#############################################################################

begcom = re.compile(r'^[^%]*\\begin\{comment\}')
endcom = re.compile(r'^[^%]*\\end\{comment\}')
isplot = re.compile(r'^[^%]*\\(?:side|full)?plot\*?\s*(?:\[[\!htbp]+\])?' \
                    '\{([^\}]+)')
ismplot = re.compile(r'^[^%]*\\multiplot\*?\s*(?:\[[\!htbp]+\])?' \
                     '\{[^\}]+\}\s*\{([^\}]+)')
issmplot = re.compile(r'^[^%]*\\sidemultiplot\*?\s*(?:\[[\!htbp]+\])?' \
                     '\{[^\}]+\}\s*\{([^\}]+)')
isfig  = re.compile(r'^[^%]*\\includegraphics\s*(\[[^\]]*\])?\{([^\}]+)')
isanim = re.compile(r'^[^%]*\\animategraphics\s*(\[[^\]]*\])?\{([0-9]+)\}\{([^\}]+)')
isbib = re.compile(r'\\bibliography\s*\{([^\}]+)')
linput = re.compile(r'[^%]\\(?:lst)?input(?:listing\[[^\]]+\])?\s*\{([^\}]+)')
chdir = re.compile(r'[^%]*\\inputdir\s*\{([^\}]+)')
subdir = re.compile(r'\\setfigdir{([^\}]+)')
beamer = re.compile(r'\\documentclass[^\{]*\{beamer\}')
hastoc =  re.compile(r'\\tableofcontents')
figure = re.compile(r'\\contentsline \{figure\}\{\\numberline \{([^\}]+)')
subfigure = re.compile(r'\\contentsline \{subfigure\}\{\\numberline \{\(([\w])')
logfigure = re.compile(r'\s*\<use ([^\>]+)')
suffix = re.compile('\.[^\.]+$')
cwpslides = re.compile(r'\\documentclass[^\{]*\{cwpslides\}')
makeind = re.compile(r'\\index')

#############################################################################
# CUSTOM SCANNERS
#############################################################################

plotoption = {}
geomanuscript = 0
slides = 0

def latexscan(node,env,path):
    'Scan LaTeX file for extra dependencies'
    global plotoption, geomanuscript, slides
    lclass = env.get('lclass','geophysics')
    options = env.get('options')
    if not options:
        geomanuscript = 0
    else:
        geomanuscript = lclass == 'geophysics' and \
            options.rfind('manuscript') >= 0
    slides = lclass == 'beamer' or lclass == 'cwpslides'

    top = str(node)
    if top[-4:] != '.ltx':
        return []
    contents = node.get_contents()
    inputs = list(filter(os.path.isfile,
                    [x+('.tex','')[os.path.isfile(x)] for x in linput.findall(contents.decode('utf-8'))]))
    inputs.append(top[:-4]+'.tex')
    resdir = env.get('resdir','Fig')
    inputdir = env.get('inputdir','.')
    plots = []
    for file in inputs:
        try:
            inp = open(file,'r')
        except:
            return []
        comment = 0
        for line in inp.readlines():
            if comment:
                if endcom.search(line):
                    comment = 0
                continue
            if begcom.search(line):
                comment = 1
                continue

            dir  = chdir.search(line)
            if dir:
                inputdir = dir.group(1)
            dir = subdir.search(line)
            if dir:
                resdir = dir.group(1)
            resdir2 = os.path.join(inputdir,resdir)

            check = isplot.search(line)
            if check:
                 plot = check.group(1)
                 plot = plot.replace('\_','_')
                 plots.append(os.path.join(resdir2,plot + ressuffix))
                 if re.search('angle=90',line):
                      plotoption[plot+pssuffix] = '-flip r90'

            check = ismplot.search(line)
            if check:
                 mplot = check.group(1)
                 mplot = mplot.replace('\_','_')
                 for plot in mplot.split(','):
                     plots.append(os.path.join(resdir2,plot + ressuffix))
                     if re.search('angle=90',line):
                         plotoption[plot+pssuffix] = '-flip r90'

            check = issmplot.search(line)
            if check:
                 smplot = check.group(1)
                 smplot = smplot.replace('\_','_')
                 for plot in smplot.split(smplot,','):
                     plots.append(os.path.join(resdir2,plot + ressuffix))
                     if re.search('angle=90',line):
                         plotoption[plot+pssuffix] = '-flip r90'

            check = isfig.search(line)
            if check:
                 plot = check.group(2)
                 if plot[-len(ressuffix):] != ressuffix:
                     plot = plot + ressuffix
                 plots.append(plot)

            check = isanim.search(line)
            if check:
                 plot = check.group(3)
                 # Only make VPL->PDF target if the file name does
                 # not end with '_', otherwise, assume that \animategraphics
                 # is going to assemble the animation from individual files
                 if plot[len(plot)-1] != '_':
                     if plot[-len(ressuffix):] != ressuffix:
                         plotoption[plot+pssuffix] = ' cropshift=y'
                         plot = plot + ressuffix
                     else:
                         plotoption[plot[-len(ressuffix):]+pssuffix] = ' cropshift=y'
                     plots.append(plot)

        inp.close()
    bibs = []
    for bib in isbib.findall(contents.decode('utf-8')):
        for file in bib.split(','):
            file = file+'.bib'
            if os.path.isfile(file):
                bibs.append(file)
    check = plots + inputs + bibs
    return check

LaTeXS = Scanner(name='LaTeX',function=latexscan,skeys=['.ltx'])

#############################################################################
# CUSTOM BUILDERS
#############################################################################

def latify(target=None,source=None,env=None):
    "Add header and footer to make a valid LaTeX file"
    tex = open(str(source[0]),'r')
    ltx = open(str(target[0]),'w')
    lclass = env.get('lclass','geophysics')
    if lclass == 'segabs':
        size = '11pt'
    else:
        size = '12pt'
    options = env.get('options',size)
    if not options:
        options = size

    ltx.write('%% This file is automatically generated. Do not edit!\n')
    ltx.write('\\documentclass[%s]{%s}\n\n' % (options,lclass))
    use = env.get('use')
    resdir = env.get('resdir','Fig')
    include = env.get('include')
    if 'endfloat' in options.split(','):
        notendfloat=False
    else:
        notendfloat=True
    if use:
         if type(use) is not list:
              use = use.split(',')
         for package in use:
              options = re.match(r'(\[[^\]]*\])\s*(\S+)',package)
              if options:
                   ltx.write('\\usepackage%s{%s}\n' % options.groups())
              else:
                   ltx.write('\\usepackage{%s}\n' % package)
         ltx.write('\n')
    if include:
        ltx.write(include+'\n\n')
    if lclass in ('segabs','georeport'):
        ltx.write('\\setfigdir{%s}\n\n' % resdir)
    if lclass == 'geophysics':
        if notendfloat:
            ltx.write('\\setfigdir{%s}\n\n' % resdir)
    ltx.write('\\begin{document}\n')
    for line in tex.readlines():
        ltx.write(line)
    ltx.write('\\end{document}\n')
    ltx.close()
    return 0

def sage_emit(target=None, source=None, env=None):
    sage = str(source[0])
    target.append(sage+'.py')
    return target, source

def latex_emit(target=None, source=None, env=None):
    tex = str(source[0])
    stem = suffix.sub('',tex)
    target.append(stem+'.aux')
    target.append(stem+'.log')
    contents = source[0].get_contents()
    if isbib.search(contents.decode('utf-8')):
        target.append(stem+'.bbl')
        target.append(stem+'.blg')
    if hastoc.search(contents.decode('utf-8')):
        target.append(stem+'.toc')
    if beamer.search(contents.decode('utf-8')):
        target.append(stem+'.nav')
        target.append(stem+'.out')
        target.append(stem+'.snm')
    if cwpslides.search(contents.decode('utf-8')):
        target.append(stem+'.nav')
        target.append(stem+'.out')
        target.append(stem+'.snm')
        target.append(stem+'.toc')
    if makeind.search(contents.decode('utf-8')):
        target.append(stem+'.idx')
    return target, source

def latex2dvi(target=None,source=None,env=None):
    "Convert LaTeX to DVI/PDF"
    tex = str(source[0])
    dvi = str(target[0])
    stem = suffix.sub('',dvi)
    if not latex:
        print('\n\tLaTeX is missing. ' \
            'Please install a TeX package (teTeX or TeX Live)\n')
        return 1
    run = ' '.join([latex,tex])
    # First latex run
    if os.system(run):
        return 1
    # Check if bibtex is needed
    aux = open(stem + '.aux',"r")
    for line in aux.readlines():
        if re.search("bibdata",line):
            if not bibtex:
                print('\n\tBibTeX is missing.')
                return 1
            os.system(' '.join([bibtex,stem]))
            os.system(run)
            os.system(run)
            break
        elif re.search("beamer@",line):
            os.system(run)
            os.system(run)
            break
    aux.close()
    # Check if makeindex is needed
    idx = stem + '.idx'
    if os.path.isfile(idx):
        if not makeindex:
            print('\n\tMakeIndex is missing.')
            return 1
        os.system(' '.join([makeindex,idx]))
        os.system(run)
    # Check if rerun is needed
    for i in range(3): # repeat 3 times at most
        done = 1
        if sys.version_info[0] >= 3:
            log = open(stem + '.log',"r",encoding='latin1')
        else:
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

loffigs = {}

def listoffigs(target=None,source=None,env=None):
    "Copy figures"
    global loffigs

    pdf = str(source[0])
    stem = suffix.sub('',pdf)

    try:
        lof = open(stem+'.lof')
        log = open(stem+'.log')
    except:
        return target, source

    figs = []
    for line in lof.readlines():
        fig = figure.match(line)
        if fig:
            figs.append(fig.group(1))
        else:
            subfig = subfigure.match(line)
            if subfig:
                last = figs.pop()
                if re.search('[a-z]$',last):
                    figs.append(last)
                    figs.append(last[:-1]+subfig.group(1))
                else:
                    figs.append(last+subfig.group(1))
    lof.close()

    for line in log.readlines():
        fil = logfigure.match(line)
        if fil and figs:
            fig = figs.pop(0)
            for ext in ('eps','pdf'):
                src = suffix.sub('.'+ext,fil.group(1))
                dst = '%s-Fig%s.%s' % (stem,fig,ext)

                loffigs[dst]=src
                target.append(dst)
                env.Precious(dst)
    log.close()

    return target, source

def copyfigs(target=None,source=None,env=None):
     "Copy figures"
     for fig in target[1:]:
          dst = str(fig)
          src = loffigs[dst]
          try:
               shutil.copy(src,dst)
          except:
               sys.stderr.write('Cannot copy %s\n' % src)
     return 0

def latex2mediawiki(target=None,source=None,env=None):
    "Convert LaTeX to MediaWiki"
    texfile = str(source[0])
    tex = open(texfile,"r")
    bblfile = re.sub('\.[^\.]+$','.bbl',texfile)
    try:
        bbl = open(bblfile,"r")
        rsf.latex2wiki.parse_bbl(bbl)
        bbl.close()
    except:
        pass
    wiki = open(str(target[0]),"w")
    rsf.latex2wiki.convert(tex,wiki)
    wiki.close()
    tex.close()
    return 0

colorfigs = []
hiresfigs = []

def pstexpen(target=None,source=None,env=None):
    "Convert vplot to EPS"
    global colorfigs, geomanuscript, plotoption

    vpl = str(source[0])
    eps = str(target[0])
    ploption = plotoption.get(eps,'')

    if vpl[-len(pssuffix):]==pssuffix:
        try:
            vplfile = open(vpl,'r')
            epsfile = open(eps,'w')
            epsfile.write(vplfile.read())
            vplfile.close()
            epsfile.close()
        except:
            sys.stderr.write('EPS write failed\n')
            return 1
    else:
        try:
            import rsf.vpconvert as vpconvert
            options = 'color=n fat=1 fatmult=1.5 invras=y'
            name = os.path.splitext(os.path.basename(eps))[0]
            if 'ALL' in colorfigs or name in colorfigs:
                options += ' color=y'
            if geomanuscript:
                options += ' serifs=n'
            elif slides:
                options += ' fat=2 txscale=1.25'
            if ploption:
                options += ploption
            vpconvert.convert(vpl,eps,'eps',None,options)
        except:
            sys.stderr.write('vpconvert failed\n')
            return 1
    return 0

_KEYWORD = token.NT_OFFSET + 1
_TEXT    = token.NT_OFFSET + 2

_colors = {
     token.NUMBER:       '#0080C0',
     token.OP:           '#0000C0',
     token.STRING:       '#004080',
     tokenize.COMMENT:   '#008000',
     token.NAME:         '#000000',
     token.ERRORTOKEN:   '#FF8080',
     _KEYWORD:           '#C00000',
     _TEXT:              '#000000',
     'Fetch':            '#0000C0',
     'Flow':             '#0000C0',
     'Plot':             '#0000C0',
     'Result':           '#C00000'
     }

_styles = {
     token.NUMBER:       'number',
     token.OP:           'op',
     token.STRING:       'string',
     tokenize.COMMENT:   'comment',
     token.NAME:         'name',
     token.ERRORTOKEN:   'error',
     _KEYWORD:           'keyword',
     _TEXT:              'text',
     'Fetch':            'fetch',
     'Flow':             'flow',
     'Plot':             'plot',
     'Result':           'result'
     }

_pos = 0


def _proglink(name):
    link = '<a href="http://www.ahay.org/RSF/%s.html">%s</a>' % (rsf.doc.progs[name].name, name)
    return link

dataserver = os.environ.get('RSF_DATASERVER','http://www.ahay.org')

def _datalink(name):
    global dataserver
    if name[:4] == 'http' or name[:3] == 'ftp':
        link = '<a href="%s">%s</a>' % (name,name)
    else:
        link = '<a href="%s/%s">%s</a>' % (dataserver,name,name)
    return link

def colorize(target=None,source=None,env=None):
     "Colorize python source"
     py = str(source[0])
     html = str(target[0])

     src = open(py,'r').read()
     raw = src.expandtabs().strip()

     out = open(html,'w')
     out.write('''
     <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN">
     <html>
     <head>
     <title>%s</title>
     <style type="text/css">
     div.progs {
     background-color: #DCE3C4;
     border: thin solid black;
     padding: 1em;
     margin-left: 2em;
     margin-right: 2em; }
     div.dsets {
     background-color: #E3C4DC;
     border: thin solid black;
     padding: 1em;
     margin-left: 2em;
     margin-right: 2em; }
     div.scons {
     background-color: #FFF8ED;
     border: thin solid black;
     padding: 1em;
     margin-left: 2em;
     margin-right: 2em; }
     ''' % py)
     for style in list(_styles.keys()):
          out.write('.%s { color: %s; }\n' % (_styles[style],_colors[style]))
     out.write('''</style>
     </head>
     <body>
     <div>
     <a href="paper_html/paper.html"><img width="32" height="32"
     align="bottom" border="0" alt="up" src="paper_html/icons/up.%s"></a>
     <a href="paper.pdf"><img src="paper_html/icons/pdf.%s" alt="[pdf]"
     width="32" height="32" border="0"></a>
     </div>
     <div class="scons">
     <table><tr><td>
     ''' % (itype,itype))

     # store line offsets in self.lines
     lines = [0, 0]
     _pos = 0
     while 1:
          _pos = raw.find('\n', _pos) + 1
          if not _pos: break
          lines.append(_pos)
     lines.append(len(raw))


     # parse the source and write it
     _pos = 0
     text = io.StringIO(unicode(raw))
     out.write('<pre><font face="Lucida,Courier New">')

     def call(toktype, toktext, xxx_todo_changeme, xxx_todo_changeme1, line):
          (srow,scol) = xxx_todo_changeme
          (erow,ecol) = xxx_todo_changeme1
          global _pos

          # calculate new positions
          oldpos = _pos
          newpos = lines[srow] + scol
          _pos = newpos + len(toktext)

          # handle newlines
          if toktype in [token.NEWLINE, tokenize.NL]:
               out.write("\n")
               return

          # send the original whitespace, if needed
          if newpos > oldpos:
               out.write(raw[oldpos:newpos])

          # skip indenting tokens
          if toktype in [token.INDENT, token.DEDENT]:
               _pos = newpos
               return

          # map token type to a color group
          if token.LPAR <= toktype and toktype <= token.OP:
               toktype = token.OP
          elif toktype == token.NAME and keyword.iskeyword(toktext):
               toktype = _KEYWORD
          elif toktype == token.NAME and toktext in list(_colors.keys()):
               toktype = toktext

          style = _styles.get(toktype, _styles[_TEXT])

          # send text
          out.write('<span class="%s">' % style)
          out.write(cgi.escape(toktext))
          out.write('</span>')

     try:
          tokenize.tokenize(text.readline, call)
     except tokenize.TokenError as ex:
          msg = ex[0]
          line = ex[1][0]
          out.write("<h3>ERROR: %s</h3>%s\n" % (msg, raw[lines[line]:]))
          return 1

     out.write('</font></pre></table>')

     info = str(source[1])

     if os.path.isfile(info):
         sout = open(info)
         progs = sout.read()
         sout.close()

         eval(progs, locals())

         if uses:
             out.write('</div><p><div class="progs">')
             out.write(rsf.doc.multicolumn(uses,_proglink))

         if data:
             while 'PRIVATE' in data:
                 data.remove('PRIVATE')
             while 'LOCAL' in data:
                 data.remove('LOCAL')
             out.write('</div><p><div class="dsets">')
             out.write(rsf.doc.multicolumn(data,_datalink))

     out.write('''
     </div>
     </body>
     </html>
     ''')

     return 0

def eps2png(target=None,source=None,env=None):
    global plotoption

    png = str(target[0])
    eps = str(source[0])
    option = plotoption.get(os.path.basename(eps),'')
    command =  'PAPERSIZE=ledger %s %s -out %s' \
              + ' -type %s -interlaced -antialias -crop a %s'
    if not pstoimg:
        print('\n\t"pstoimg" is missing. ' \
            'Please install latex2html.\n')
        return 1
    command = command % (pstoimg,eps,png,itype,option)
    print(command)
    os.system(command)
    return 0

def eps2pdf(target=None,source=None,env=None):
    global hiresfigs

    pdf = str(target[0])
    eps = str(source[0])

    gs_options = os.environ.get('GS_OPTIONS','')
    if not gs_options:
        name = os.path.splitext(os.path.basename(eps))[0]
        if name in hiresfigs:
            gs_options = '-dAutoFilterColorImages=false -dColorImageFilter=/LZWEncode ' + \
                         '-dAutoFilterGrayImages=false  -dGrayImageFilter=/LZWEncode'
    if gs_options:
        gs_options = "GS_OPTIONS='%s'" % gs_options

    command = "LD_LIBRARY_PATH=%s %s %s --hires %s" % \
              (os.environ.get('LD_LIBRARY_PATH',''),gs_options,epstopdf,eps)

    print(command)
    os.system(command)
    return 0

def dummy(target=None,source=None,env=None):
     tex = open(str(target[0]),'w')
     tex.write('%% This file is automatically generated. Do not edit!\n')

     user = os.getuid()
     name = pwd.getpwuid(user)[4]

     tex.write('\\author{%s}\n' % name)
     tex.write('\\title{Dummy paper}\n\n\maketitle\n')
     dirold = ''
     for src in source:
         fig = str(src)
         plt = os.path.splitext(os.path.basename(fig))[0]
         plt2 = plt.replace('_','\_')
         fdir = os.path.split(os.path.split(fig)[0])[0]
         if fdir != dirold:
             tex.write('\n\\section{%s}\n' % fdir)
             if fdir:
                 tex.write('\\inputdir{%s}\n\n' % fdir)
             else:
                 tex.write('\\inputdir{.}\n\n')
             dirold = fdir
         tex.write('\\plot{%s}{width=\\textwidth}{%s/%s} ' % (plt,fdir,plt2))
         tex.write('\\clearpage\n')
     tex.close()
     return 0

def pylab(target=None,source=None,env=None):
    global epstopdf
    pycomm = open(str(source[0]),'r').read()
    exec(pycomm, locals())
    os.system('%s junk_py.eps -o=%s' % (epstopdf,target[0]))
    os.unlink('junk_py.eps')
    return 0

Latify = Builder(action = Action(latify,
                                 varlist=['lclass','options','use',
                                          'include','resdir']),
                 src_suffix='.tex',suffix='.ltx')
Pdf = Builder(action=Action(latex2dvi,varlist=['latex','lclass',
                                               'options','resdir']),
              source_scanner=LaTeXS,emitter=latex_emit,
              src_suffix='.ltx',suffix='.pdf')
Wiki = Builder(action=Action(latex2mediawiki),src_suffix='.ltx',suffix='.wiki')
Figs = Builder(action=Action(copyfigs),
               src_suffix='.pdf',emitter=listoffigs)

if pdfread:
    Read = Builder(action = pdfread + " $SOURCES",
                   src_suffix='.pdf',suffix='.read')
    Print = Builder(action =
                    'cat $SOURCES | %s -toPostScript | lpr' % pdfread,
                    src_suffix='.pdf',suffix='.print')

Build = Builder(action = Action(pstexpen),
                src_suffix=vpsuffix,suffix=pssuffix)

if epstopdf:
    PDFBuild = Builder(action = Action(eps2pdf),
		       src_suffix=pssuffix,suffix='.pdf')

if fig2dev:
    XFig = Builder(action = fig2dev + ' -L pdf -p dummy $SOURCES $TARGETS',
                   suffix='.pdf',src_suffix='.fig')

PNGBuild = Builder(action = Action(eps2png),
                   suffix='.'+itype,src_suffix=pssuffix)

if pdftops:
    PSBuild = Builder(action = pdftops + ' -eps $SOURCE $TARGET',
                      suffix=pssuffix,src_suffix='.pdf')
elif acroread and ps2eps:
    gs = WhereIs('gs')
    if gs:
        PSBuild = Builder(action = '%s -toPostScript -size ledger -pairs $SOURCE'
                          ' junk.ps && %s -q -sDEVICE=epswrite -sOutputFile=$TARGET'
                          ' -r600 -dNOPAUSE -dBATCH junk.ps && rm junk.ps' % \
                          (acroread,gs),
                          suffix=pssuffix,src_suffix='.pdf')
    else:
        PSBuild = Builder(action = '%s -toPostScript -size ledger -pairs $SOURCE'
                          ' junk.ps && %s junk.ps $TARGET && rm junk.ps' % \
                          (acroread,ps2eps),
                          suffix=pssuffix,src_suffix='.pdf')
elif pdf2ps:
    PSBuild = Builder(action = pdf2ps + ' $SOURCE $TARGET',
                      suffix=pssuffix,src_suffix='.pdf')

if latex2html:
    l2hdir = os.environ.get('LATEX2HTML','')
    inputs = os.environ.get('TEXINPUTS','')

    if not l2hdir:
        init = ''
        proc = os.popen(WhereIs('kpsewhich') + ' -show-path ls-R 2>/dev/null')
        raw  = proc.read()
        if raw:
            for dir in raw.split(':'):
                ldir = os.path.join(dir,'latex2html')
                if os.path.isdir(ldir):
                    l2hdir = ldir
                    break
    if l2hdir:
        init = '-init_file ' + os.path.join(l2hdir,'.latex2html-init')
        css0 = os.path.join(l2hdir,'style.css')
        icons0 = os.path.join(l2hdir,'icons')

    HTML = Builder(action = 'TEXINPUTS=%s LATEX2HTMLSTYLES=%s/perl %s '
                   '-debug $SOURCE -dir $TARGET.dir -image_type %s %s' %
                   (inputs,l2hdir,latex2html,itype,init),src_suffix='.ltx')

if sage:
    Sage = Builder(action = sage + ' $SOURCE && mv junk_sage.pdf $TARGET',
                   suffix='.pdf',src_suffix='.sage',emitter=sage_emit)

if epstopdf:
    if mathematica:
        Math = Builder(action = '%s -batchoutput '
                       '< $SOURCE >&2  > /dev/null && '
                       '%s --hires junk_ma.eps -o=$TARGET && rm junk_ma.eps' %
                       (mathematica,epstopdf),
                       suffix='.pdf',src_suffix='.ma')
    if matlab:
        matlabpath = os.environ.get('MATLABPATH')
        if matlabpath:
            matlabpath = ':'.join([matlabpath,'Matlab'])
        else:
            matlabpath = 'Matlab'
        Matlab = Builder(action = 'MATLABPATH=%s DISPLAY=" " nohup %s -nodesktop '
                         '< $SOURCE >&2  > /dev/null && '
                         '%s junk_ml.eps -o=$TARGET && rm junk_ml.eps' %
                         (matlabpath,matlab,epstopdf),
                         suffix='.pdf',src_suffix='.ml')
    if gnuplot:
        Gnuplot = Builder(action = '%s $SOURCE > junk_gp.eps && '
                          '%s junk_gp.eps -o=$TARGET && rm junk_gp.eps' %
                          (gnuplot,epstopdf),
                          suffix='.pdf',src_suffix='.gp')

    if haspylab:
        Pylab = Builder(action = Action(pylab),
                        suffix='.pdf',src_suffix='.py')

Color = Builder(action = Action(colorize),suffix='.html')

class TeXPaper(Environment):
    def __init__(self,**kw):
        kw.update({'tools':[]})
        Environment.__init__(*(self,), **kw)
        rsf.conf.set_options(self)
#        sourceforge = 'http://sourceforge.net/p/rsf/code/HEAD/tree/trunk'
        github = 'https://github.com/ahay/src/blob/master/'
        self.Append(ENV={'XAUTHORITY':
                         os.path.join(os.environ.get('HOME'),'.Xauthority'),
                         'DISPLAY': os.environ.get('DISPLAY'),
			 'RSF_REPOSITORY': os.environ.get('RSF_REPOSITORY',github),
			 'RSF_ENSCRIPT': WhereIs('enscript'),
                         'HOME': os.environ.get('HOME')},
                    SCANNERS=LaTeXS,
                    BUILDERS={'Latify':Latify,
                              'Pdf':Pdf,
                              'Wiki':Wiki,
                              'Build':Build,
                              'Color':Color,
                              'Figs':Figs})
        path = {'darwin': ['/sw/bin','/opt/local/bin'],
                'irix': ['/usr/freeware/bin']}
        for plat in list(path.keys()):
            if sys.platform[:len(plat)] == plat:
                for pathdir in filter(os.path.isdir,path[plat]):
                    self['ENV']['PATH'] = ':'.join([pathdir,
                                                    self['ENV']['PATH']])

        tree = rsf.path.dirtree()

        root = self.get('RSFROOT',rsf.prog.RSFROOT)
        self.docdir = os.environ.get('RSFBOOK',os.path.join(root,'share','madagascar'))
        self.figdir = os.environ.get('RSFFIGS',os.path.join(self.docdir,'figs'))

        for level in tree:
            if level:
                self.docdir = os.path.join(self.docdir,level)
        rsf.path.mkdir(self.docdir)

        datapath = rsf.path.datapath()
        self.path = os.path.dirname(datapath)
        if datapath[:2] != './':
            for level in tree[1:]:
                self.path = os.path.join(self.path,level)
        rsf.path.mkdir(self.path)
        self.path = os.path.join(self.path,os.path.basename(datapath))
        rsf.path.sconsign(self)

        if pdfread:
            self.Append(BUILDERS={'Read':Read,'Print':Print})
        if epstopdf:
            self.Append(BUILDERS={'PDFBuild':PDFBuild})
        if fig2dev:
            self.Append(BUILDERS={'XFig':XFig})
        if latex2html:
            self.Append(BUILDERS={'HTML':HTML,'PNGBuild':PNGBuild})
            self.imgs = []
        if (acroread and ps2eps) or pdf2ps:
            self.Append(BUILDERS={'PSBuild':PSBuild})
        if epstopdf:
            if mathematica:
                self.Append(BUILDERS={'Math':Math})
            if gnuplot:
                self.Append(BUILDERS={'Gnuplot':Gnuplot})
            if matlab:
                self.Append(BUILDERS={'Matlab':Matlab})
            if haspylab:
                self.Append(BUILDERS={'Pylab':Pylab})
            if sage:
                self.Append(BUILDERS={'Sage':Sage})

        self.scons = []
        self.figs = []
        self.Dir()
    def Install2(self,dir,fil):
        dir2 = rsf.path.mkdir(dir)
        self.Install(dir2,fil)
    def Dir(self,topdir='.',resdir='Fig'):
        # reproducible directories
        for info in glob.glob('%s/[a-z]*/.*proj' % topdir):
            dir = os.path.dirname(info)
            scons = os.path.join(dir,'SConstruct')

            html = dir+'.html'
            self.Color(html,[scons,info])
            self.scons.append(html)

        if self.scons:
            self.Install(self.docdir,self.scons)
        self.Alias('figinstall',self.docdir)
        # reproducible figures
        erfigs = []
        eps = {}


        # check figure repository
        vpldir = re.sub(r'.*\/((?:[^\/]+)\/(?:[^\/]+))$',
                        self.figdir+'/\\1',os.path.abspath(topdir))
        for suffix in (vpsuffix,pssuffix):
            for fig in glob.glob('%s/[a-z]*/*%s' % (vpldir,suffix)):
                eps[fig] = re.sub(r'.*\/([^\/]+)\/([^\/]+)'+suffix+'$',
                                  r'%s/\1/%s/\2%s' % (topdir,resdir,pssuffix),
                                  fig)

        # follow symbolic links
        for pdir in filter(os.path.islink,glob.glob(topdir+'/[a-z]*')):
            vpldir = re.sub(r'.*\/((?:[^\/]+)\/(?:[^\/]+)\/(?:[^\/]+))$',
                            self.figdir+'/\\1',
                            os.path.abspath(os.path.realpath(pdir)))
            for suffix in (vpsuffix,pssuffix,'.pdf'):
                for fig in glob.glob('%s/*%s' % (vpldir,suffix)):
                    eps[fig] = re.sub(r'.*\/([^\/]+)\/([^\/]+)'+suffix+'$',
                                      r'%s/%s/\2%s' % (pdir,resdir,pssuffix),
                                      fig)

        for fig in list(eps.keys()):
            ps = eps[fig]
            resdir2 = os.path.join(self.docdir,os.path.dirname(ps))
            if fig[-3:] == vpsuffix[-3:]:
                self.Build(ps,fig)
            else:
                self.InstallAs(ps,fig)
            if epstopdf:
                pdf = re.sub(pssuffix+'$','.pdf',ps)
                self.PDFBuild(pdf,ps)
                erfigs.append(pdf)
                self.Install2(resdir2,pdf)
            if latex2html and pstoimg:
                png = re.sub(pssuffix+'$','.'+itype,ps)
                self.PNGBuild(png,ps)
                self.imgs.append(png)
                self.Install2(resdir2,png)
                self.Alias('figinstall',resdir2)
        self.figs.extend(erfigs)

        # conditionally reproducible figures
        crfigs = []
        # mathematica figures:
        mths = glob.glob('%s/Math/*.ma' % topdir)
        if mths:
            for mth in mths:
                pdf = re.sub(r'([^/]+)\.ma$',
                             os.path.join(resdir,'\g<1>.pdf'),mth)
                if mathematica and epstopdf:
                    self.Math(pdf,mth)
                crfigs.append(pdf)
            mathdir = os.path.join(self.docdir,'Math')
            self.Install2(mathdir,mths)
            self.Alias('figinstall',mathdir)
        # gnuplot figures:
        gpls = glob.glob('%s/Gnuplot/*.gp' % topdir)
        if gpls:
            for gpl in gpls:
                pdf = re.sub(r'([^/]+)\.gp$',
                             os.path.join(resdir,'\g<1>.pdf'),gpl)
                if gnuplot and epstopdf:
                    self.Gnuplot(pdf,gpl)
                crfigs.append(pdf)
            gpldir = os.path.join(self.docdir,'Gnuplot')
            self.Install2(gpldir,gpls)
            self.Alias('figinstall',gpldir)
        # sage figures
        sapls = glob.glob('%s/Sage/*.sage' % topdir)
        if sapls:
            for sapl in sapls:
                pdf = re.sub(r'([^/]+)\.sage$',
                             os.path.join(resdir,'\g<1>.pdf'),sapl)
                if sage and epstopdf:
                    self.Sage(pdf,sapl)
                crfigs.append(pdf)
            sapldir = os.path.join(self.docdir,'Sage')
            self.Install2(sapldir,sapls)
            self.Alias('figinstall',sapldir)
        # matlab figures
        mtls = glob.glob('%s/Matlab/*.ml' % topdir)
        if mtls:
            for mtl in mtls:
                pdf = re.sub(r'([^/]+)\.ml$',
                             os.path.join(resdir,'\g<1>.pdf'),mtl)
                if matlab and epstopdf:
                    self.Matlab(pdf,mtl)
                crfigs.append(pdf)
            matlabdir = os.path.join(self.docdir,'Matlab')
            self.Install2(matlabdir,mtls)
            self.Alias('figinstall',matlabdir)
        # pylab figures
        pyls = glob.glob('%s/Pylab/*.py' % topdir)
        if pyls:
            for pyl in pyls:
                pdf = re.sub(r'([^/]+)\.py$',
                             os.path.join(resdir,'\g<1>.pdf'),pyl)
                if haspylab and epstopdf:
                    self.Pylab(pdf,pyl)
                crfigs.append(pdf)
            pylabdir = os.path.join(self.docdir,'Pylab')
            self.Install2(pylabdir,pyls)
            self.Alias('figinstall',pylabdir)
        # xfig figures:
        figs =  glob.glob('%s/XFig/*.fig' % topdir)
        if figs:
            for fig in figs:
                pdf = re.sub(r'([^/]+)\.fig$',
                             os.path.join(resdir,'\g<1>.pdf'),fig)
                if fig2dev:
                    self.XFig(pdf,fig)
                crfigs.append(pdf)
            resdir2 = os.path.join(self.docdir,'XFig')
            self.Install2(resdir2,figs)
            self.Alias('figinstall',resdir2)
        # tikz figures:
        figs =  glob.glob('%s/Tikz/*.tex' % topdir)
        if figs:
            for fig in figs:
                pdf = re.sub(r'([^/]+)\.tex$',
                             os.path.join(resdir,'\g<1>.pdf'),fig)
                crfigs.append(pdf)
            resdir2 = os.path.join(self.docdir,'Tikz')
            self.Install2(resdir2,figs)
            self.Alias('figinstall',resdir2)
        # non-reproducible figures
        nrfigs = crfigs + glob.glob(
            os.path.join(topdir,os.path.join(resdir,'*.pdf')))
        for pdf in nrfigs:
             if (acroread and ps2eps) or pdf2ps:
                eps = re.sub('.pdf$',pssuffix,pdf)
                self.PSBuild(eps,pdf)
                if latex2html and pstoimg:
                    png = re.sub(pssuffix+'$','.'+itype,eps)
                    self.PNGBuild(png,eps)
                    self.imgs.append(png)
                    resdir2 = os.path.join(self.docdir,os.path.dirname(png))
                    self.Install2(resdir2,[png,pdf])
                    self.Alias('figinstall',resdir2)
        self.figs.extend(nrfigs)
        self.Command('dummy.tex',self.figs,Action(dummy))
    def Paper(self,paper,source='',lclass='geophysics',scons=1,
              use=None,include=None,options=None,
              resdir='Fig',color='',hires=''):
        global colorfigs, hiresfigs
        if source == '':
            source = paper
        colorfigs.extend(color.split())
        hiresfigs.extend(hires.split())
        ltx = self.Latify(target=paper+'.ltx',source=source+'.tex',
                          use=use,lclass=lclass,options=options,
                          include=include,resdir=resdir)
        pdf = self.Pdf(target=paper,source=paper+'.ltx',
                       lclass=lclass,options=options,resdir=resdir)
        self.Figs(target=paper+'.figs',source=paper+'.pdf')
        wiki = self.Wiki(target=paper,source=[ltx,pdf])
        pdfinstall = self.Install(self.docdir,paper+'.pdf')
        self.Alias(paper+'.install',pdfinstall)
        if pdfread:
            self.Alias(paper+'.read',self.Read(paper))
            self.Alias(paper+'.print',self.Print(paper))
        if latex2html and l2hdir:
            hdir = paper+'_html'
            css  = os.path.join(hdir,paper+'.css')
            html = os.path.join(hdir,'index.html')
            icons = os.path.join(hdir,'icons')
            self.InstallAs(css,css0)
            self.Install(icons,glob.glob('%s/*.%s' % (icons0,itype)))
            self.HTML(html,paper+'.ltx')
            self.Depends(self.imgs,pdf)
            self.Depends(html,self.imgs)
            if scons:
                self.Depends(html,self.scons)
            self.Depends(html,pdf)
            self.Depends(html,css)
            self.Depends(html,icons)
            self.Alias(paper+'.html',html)
            docdir = os.path.join(self.docdir,hdir)
            dochtml = os.path.join(docdir,'index.html')
            self.Command(dochtml,html,
                         'cd $SOURCE.dir && cp -R * $TARGET.dir && cd ..')
            self.Alias(paper+'.install',dochtml)
            self.Depends(paper+'.install','figinstall')
        return pdf
    def End(self,paper='paper',source='',**kw):
        if source == '':
            source = paper
        if os.path.isfile(source+'.tex'):
            self.Paper(*(paper,source), **kw)
            self.Alias('pdf',paper+'.pdf')
            self.Alias('wiki',paper+'.wiki')
            self.Alias('read',paper+'.read')
            self.Alias('print',paper+'.print')
            self.Alias('html',paper+'.html')
            self.Alias('install',paper+'.install')
            self.Alias('figs',paper+'.figs')
            self.Default('pdf')

default = TeXPaper()
def Dir(**kw):
     return default.Dir(*[], **kw)
def Paper(paper,source='',**kw):
    return default.Paper(*(paper,source), **kw)
def Command2(target,source,command):
    return default.Command(target,source,command)
def End(paper='paper',source='',**kw):
    return default.End(*(paper,source), **kw)
def Depends2(target,source):
    return default.Depends(target,source)

if __name__ == "__main__":
     import pydoc
     pydoc.help(TeXPaper)


