#!/usr/bin/env python
#
# Original idea and some code from Maxime Biais. Significant contributions by
# Marc Poulhiès. Additions and changes by Mike Meylan, Malte, Anthony Miller and
# Sergey Fomel.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import re, string, os

bdoc = None
verbatim_mode = 0
math_mode = 0
in_math = 0

el = "\n"
item = ""
refs = {}
gref = []
fullref = None
insert = ""

lang = "text"
code = ""
coderef = ""
repos = "http://sourceforge.net/p/rsf/code/HEAD/tree/trunk/"

def dummy(s):
    pass

def items(s):
    global item, el
    item = item + '*'
    el = " "

def enums(s):
    global item, el
    item = item + '#'
    el = " "

def in_list():
    global item
    if item:
        return "\n%s:" % item
    else:
        return ""

def end_list(s):
    global item, el
    item = item[:-1] # remove last char
    if not item:
        el = "\n"
	
def start_doc(s):
    global bdoc;
    bdoc = 1

def decide_el():
    global el
    return el

def decide_math():
    global math_mode, verbatim_mode
    if verbatim_mode:
        math_mode = 0
        return "&#36;"
    elif math_mode:
        return "<math>"
    else:
        return "</math>"

def line_math():
    global math_mode, verbatim_mode, in_math
    if verbatim_mode:
        return r"&#36;\1&#36;"
    else:
        in_math = 1
        return r"<math>\1</math>"
    
def start_verbatim(s):
    global verbatim_mode
    verbatim_mode = 1

def end_verbatim(s):
    global verbatim_mode
    verbatim_mode = 0

def toggle_math(s):
    global math_mode, in_math
    math_mode = 1 - math_mode
    in_math = 1

def cite(s):
    global refs, gref, fullref
    fullref = s.group(1)
    keys = s.group(2)
    gref = []
    for key in string.split(keys,','):
        lref = refs.get(key)
        if lref:
            gref.append(lref)

def refer():
    global gref, fullref
    list = []
    name = ''
    for lref in gref:
        if fullref:
            list.append('%s<ref>%s</ref>' % lref)
        else:
            (name,year) = string.split(lref[0],', ')
            list.append('%s<ref>%s</ref>' % (year,lref[1]))
    lrefs = '(%s)' % string.join(list,';')
    if fullref:
        return lrefs
    else:
        return '%s %s' % (name,lrefs)

def input_file(s):
    global insert
    name = s.group(1)
    try:
        inp = open(name+'.wiki','r')
        insert = string.join(inp.readlines(),'')
        inp.clos()
    except:
        inp = '[name]'

def insert_file():
    global insert
    return insert

def blank():
    global verbatim_mode, math_mode
    if verbatim_mode or math_mode:
        return r"\1"
    else:
        return ""

def brace():
    global verbatim_mode, math_mode
    if verbatim_mode or math_mode:
        return r"\1"
    else:
        return r"\2"

def lstset(s):
    global lang
    lang = s.group(1)
    if lang == "c++":
        lang = "cpp"

def putcode():
    global lang, code
    return "<syntaxhighlight lang=\"%s\">\n%s</syntaxhighlight>\n" % (lang,code)

def getcode(s):
    global code
    options = s.group(1)
    name = string.replace(s.group(2),'\\RSF',os.environ.get('RSFSRC'))
    line = {'first':1,'last':9999}
    for mark in line.keys():
        if options:
            match = re.search('%sline=(\d+)' % mark,options)
            if match:
                line[mark] = int(match.group(1))
    fd = open(name,"r")
    code = ''
    num = 1
    for each in fd.readlines():
        if num > line['last']:
            break
        if num >= line['first']:
            code = code + each
        num = num + 1
    fd.close()

def getmode(s):
    global code, coderef, lang
    lang = 'c'
    coderef = os.path.join(s.group(4),s.group(1))+'.'+lang
    name = os.path.join(os.environ.get('RSFSRC'),coderef)
    first = int(s.group(2))
    last = int(s.group(3))
    fd = open(name,"r")
    code = ''
    num = 1
    for each in fd.readlines():
        if num > last:
            break
        if num >= first:
            code = code + each
        num = num + 1
    fd.close()

braces = r"\{((?:[^\}\{]*)(?:[^\{\}]*\{[^\}]*\}[^\{\}]*)*)\}"
# matches the contents of {} allowing a single level of nesting

tr_list2 = [
    (r"(``|'')", (lambda : '"'), dummy),
    (r"\\input{(\w+)}",insert_file,input_file),
    (r"^(\s+)", blank, dummy),
    (r"\\footnote"+braces,(lambda : r"<ref>\1</ref>"),dummy),
    (r"\\bibliography{[^}]+}",
     (lambda : r"==References==\n<references/>"),
     dummy),
    (r"\\bibliographystyle{[^}]+}",None,dummy),
    (r"\\begin\{abstract}", None, dummy),
    (r"\\begin\{article}", None, dummy),
    (r"\\end\{abstract}", None, dummy),
    (r"\\end\{article}", None, dummy),
    (r"\\end\{document}", None, dummy),
    (r"\\begin\{document}", None, start_doc),
    (r"\\inputdir{(\w+)}",None,dummy),
    (r"\\(?:side)?plot{([^}]+)}{[^}]*}{(.*)}",
     (lambda : r"[[Image:\1.png|frame|center|\2]]"),dummy),
    (r"\\label{(.*?)}", (lambda :r" (\1)"), dummy),
    (r"\\ref{(.*?)}", (lambda :r"(\1)"), dummy),
    (r"\\emph{(.*?)}", (lambda :r"''\1'' "), dummy),
    (r"\\textit{(.*?)}", (lambda :r"''\1'' "), dummy),
    (r"\\texttt{(.*?)}", (lambda : r"<tt>\1</tt>"), dummy),
    (r"\\text{(.*?)}", (lambda : r"=\1= "), dummy),
    (r"\\textbf{(.*?)}", (lambda : r"'''\1''' "), dummy),
    (r"\\verb(.)(.+)\1", (lambda : r'<font color="#cd4b19">\2</font>'), dummy),
    (r"\\begin{verbatim}", (lambda : "<pre>"), start_verbatim),
    (r"\\end{verbatim}", (lambda : "</pre>"), end_verbatim),
    (r"\\begin{comment}", (lambda : "<!-- "), dummy),
    (r"\\end{comment}", (lambda : " -->"), dummy),
    (r"\\begin{itemize}", (lambda : "\n"), items),
    (r"\\end{itemize}", in_list, end_list),
    (r"\\begin{enumerate}", (lambda : "\n"), enums),
    (r"\\end{enumerate}", in_list, end_list),
    (r"\\item (.*?)", (lambda :   r"\n" + item + r"\1"), dummy),
    (r"\\lstset{[^\}]*language=([\w\+]+)[^\}]*}",None,lstset),
    (r"\\lstinputlisting(?:\[([^\]]*)\])?{([^\}]+)}", putcode, getcode),
    (r"\\moddex{([^\}]+)}{(?:[^\}]+)}{([^\}]+)}{([^\}]+)}{([^\}]+)}",
     putcode,getmode),
    (r"\\begin{equation[*]*}", (lambda :"<center><math>"), toggle_math),
    (r"\\end{equation[*]*}", (lambda :"</math></center>"), toggle_math),
    (r"\\\[", (lambda :"<center><math>"), toggle_math),
    (r"\\dfrac", (lambda :r"\\frac"), dummy),
    (r"\\\]", (lambda :"</math></center>"), toggle_math),
    (r"\\begin{eqnarray[*]?}", (lambda :r"<center><math>\\begin{matrix}"), toggle_math),
    (r"\\begin{array[*]?}", (lambda :r"\\begin{matrix}"), toggle_math),
    (r"\\end{eqnarray[*]?}", (lambda :r"\\end{matrix}</math></center>"), toggle_math),
    (r"\\end{array[*]?}", (lambda :r"\\end{matrix}"), toggle_math),
    #	(r"(\\begin{.*?})", decide_math_replace, dummy),
    #	(r"(\\end{.*?})",decide_math_replace, dummy),
    (r"\\cite(\[\])?\{([^\}]+)\}",refer,cite),
    (r"~\\ref{([^}]*)}",(lambda : r" ---\1---"),dummy),
    (r"\\subsubsection{(.*?)}", (lambda : r"====\1===="), dummy),
    (r"\\subsection{(.*?)}", (lambda : r"===\1==="), dummy),
    (r"\\section{(.*?)}", (lambda : r"==\1=="), dummy),
    (r"\\_", (lambda :"_"), dummy),
    #	(r"\\title{(.*)}", (lambda :r"= \1 ="),dummy),
    #        (r"\\author{(.*)}", (lambda :r"\1"),dummy),
    (r"\\date{(.*)}", (lambda :r"\1"),dummy),
    (r"\\begin{quote}",(lambda : r"<blockquote>"),dummy),
    (r"\\end{quote}",(lambda : r"</blockquote>"),dummy),
    (r"\\thispagestyle{.*?}", None, dummy),
    (r"\n$", decide_el, dummy),
    #	(r"[^\\]?\{", None, dummy),
    #	(r"[^\\]?\}", None, dummy),
    (r"\\\$",(lambda : r"&#36;"),dummy),
    (r"\\\'[\{]?a[\}]?",(lambda : r"&#225;"),dummy),
    #	(r"\$(.*?)\$",(lambda :r"<math>\1</math>"),dummy),
    (r"%.*$",None, dummy),
    (r"\\r{(.*?)}", (lambda : r"\\mathrm{\1}"), dummy),
    (r"\\d ", (lambda : r"\\,\mathrm{d} "), dummy),
    (r"\\i ", (lambda : r"\\mathrm{i} "), dummy),
    (r"\\i\\", (lambda : r"\\mathrm{i}\\"), dummy),
    (r"\\e\^", (lambda : r"\\mathrm{e}^"), dummy),
    (r"\\begin{align[*]?}", (lambda :r"<center><math>\\begin{matrix}"), toggle_math),
    (r"\\end{align[*]?}", (lambda :r"\\end{matrix}</math></center>"), toggle_math),
    (r"\\begin{aligned[*]?}", None, dummy),
    (r"\\end{aligned[*]?}", None, dummy),
    (r"\\begin{subequations[*]?}", None, dummy),
    (r"\\end{subequations[*]?}", None, dummy),
    (r"\\href{([^\}]*)}{([^\}]*)}",(lambda : r"[\1 \2]"), dummy),
    (r"\\url{([^\}]*)}",(lambda : r"\1"), dummy),
    (r"\\pdfbookmark\[[^\]]*\]{[^\}]*}{[^\}]*}", None, dummy),
    # the most important thing
    (r"\\LaTeX",(lambda : "L<sup>A</sup>TEX"), dummy),
    (r"\\\\", (lambda : "\n"), dummy),
    (r"\\ ",(lambda: " "), dummy),
    (r"\$([^\$]+)\$",line_math,toggle_math),
    (r"\$",decide_math,toggle_math),
    # unknown command
    (r"(\\\w+)",blank, dummy),
    # remove braces
    (r"([{}](\w))",brace, dummy),
    (r"((\w)[{}])",brace, dummy),
    (r"^[{}]+$",None, dummy),
    ]

# precompile regular expressions
reg = map(lambda x: (re.compile(x[0]),x[1],x[2]),tr_list2)

bibitem = re.compile(r'\\bibitem\[([^\]]+)\]{([^}]+)}\s*\n(.+)$',re.DOTALL)
it_in = re.compile(r'{\\it in}')
bf = re.compile(r'{\\bf (\w+)}')
blank = re.compile(r'\n')
tilde = re.compile(r'[~]')

def parse_bbl(bbl):
    "Parse a bbl file extracting a reference dictionary"
    global refs
    for par in string.split(string.join(bbl.readlines(),''),'\n\n'):
        ref = bibitem.match(par)
        if ref:
            short = ref.group(1)
            short = tilde.sub(' ',short)

            key = ref.group(2)

            llong = ref.group(3)
            llong = it_in.sub("''in''",llong)
            llong = bf.sub("'''\\1'''",llong)
            llong = tilde.sub(' ',llong)
            llong = blank.sub('',llong)
            
            refs[key] = (short,llong)
        
def convert(in_stream,out_stream):
    "Convert LaTeX to MediaWiki"
    global reg, math_mode, in_math
    for i in in_stream.readlines():
	mystr = i

	for r in reg:
            s = r[0].search(mystr)
            if s:
                r[2](s)
            if r[1]:
                mysub = r[1]()
            else:
                mysub = ""                
            mystr = r[0].sub(mysub, mystr)
            if in_math:
                in_math = 0
                break

	if bdoc:
            out_stream.write(mystr)

if __name__ == "__main__":
    import sys, glob

    for bblfile in glob.glob('*.bbl'):
        try:
            bbl = open(bblfile,'r')
            parse_bbl(bbl)
            bbl.close()
        except:
            pass
    
    convert(sys.stdin,sys.stdout)
