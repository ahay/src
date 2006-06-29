#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
#
# Modified by Sergey Fomel <sergey.fomel@beg.utexas.edu>
#
# Mike Meylan is hacking this to include more latex commands
# and Malte has added a few things that get rid of his personal latex commands.
#
#This code has been modified by Anthony Miller for handling of inline mathematics and
#	more sophisticated documents. 
#
#Original idea from : 
#       Maxime Biais <maxime@biais.org>
#     but has been nearly all rewritten since...
# A good fraction of this code was written by
#Marc Poulhiès <marc.poulhies@epfl.ch>
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
#
# $Id: latex2twiki.py,v 1.2 2005/07/27 12:40:53 poulhies Exp $


import sys, re

bullet_level=0
enum_level = 0
bdoc = None

verbatim_mode = 0
math_mode = 0
eqnarry_mode = 0
el = "\n"

def dummy():
    pass

def inc_bullet():
	global bullet_level, el
	bullet_level += 1
        el = " "

def dec_bullet():
	global bullet_level, el
	bullet_level -= 1
        el = "\n"

def inc_enum():
	global enum_level, el
	enum_level += 1
        el = " "
	
def dec_enum():
	global enum_level, el
	enum_level -= 1
        el = "\n"
	
def start_doc():
	global bdoc;
	bdoc = 1

def decide_el():
    global el
    return el

def decide_math_replace():
	global math_mode
	if math_mode:
		return r"\1"
	else:
		return " "

def decide_math():
    global math_mode
    if verbatim_mode:
        return "&#36;"
    elif math_mode:
        return "<math>"
    else:
        return "</math>"
		
def start_verbatim():
	global verbatim_mode
	verbatim_mode = 1

def end_verbatim():
	global verbatim_mode
	verbatim_mode = 0

def start_comment():
	global bdoc
	bdoc = 0

def end_comment():
	global bdoc
	bdoc = 1

def start_eqnarry():
	global eqnarry_mode
	eqarry_mode = 1

def end_eqnarry():
	global eqnarry_mode
	eqnarry_mode = 0

def toggle_math():
	global math_mode
	math_mode = 1 - math_mode

tr_list2 = [
        (r"^\s+", None, dummy),
	(r"\\footnotesize", None, dummy),
	(r"\\begin\{abstract}", None, dummy),
	(r"\\begin\{article}", None, dummy),
	(r"\\end\{abstract}", None, dummy),
	(r"\\end\{article}", None, dummy),
	(r"\\end\{document}", None, dummy),
	(r"\\protect", None, dummy),
	(r"\\small", None, dummy),
	(r"\\func", None, dummy),
	(r"\\begin\{document}", None, start_doc),
	(r"\\cite{(.*?)}", (lambda :r"[[\1]]"), dummy),
	(r"\\label{(.*?)}", (lambda :r" (\1)"), dummy),
	(r"\\ref{(.*?)}", (lambda :r"(\1)"), dummy),
	(r"\\citep{(.*?)}", (lambda :r"[[\1]]"), dummy),
	(r"\\citet{(.*?)}", (lambda :r"[[\1]]"), dummy),
        (r"(``|'')", (lambda : '"'), dummy),
	(r"\\emph{(.*?)}", (lambda :r"''\1'' "), dummy),
	(r"\\textit{(.*?)}", (lambda :r"''\1'' "), dummy),
	(r"\\texttt{(.*?)}", (lambda : r"<tt>\1</tt>"), dummy),
	(r"\\text{(.*?)}", (lambda : r"=\1= "), dummy),
	(r"\\textbf{(.*?)}", (lambda : r"'''\1''' "), dummy),
        (r"\\verb(.)(.+)\1", (lambda : r"\2"), dummy),
	(r"\\begin{verbatim}", (lambda : "<pre>"), start_verbatim),
	(r"\\end{verbatim}", (lambda : "</pre>"), end_verbatim),
        (r"\\begin{comment}", (lambda : "<!-- "), dummy),
        (r"\\end{comment}", (lambda : " -->"), dummy),
	(r"\\begin{itemize}", (lambda : "\n"), inc_bullet),
	(r"\\end{itemize}", None, dec_bullet),
	(r"\\begin{enumerate}", (lambda : "\n"), inc_enum),
	(r"\\end{enumerate}", None, dec_enum),
	(r"\\item (.*?)", (lambda :   r"\n" + ("#" * enum_level) + \
                           ("*" * bullet_level) + r"\1"), dummy),
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
	(r"~\\ref{([^}]*)}",(lambda : r" ---\1---"),dummy),
	(r"\\subsubsection{(.*?)}", (lambda : r"====\1===="), dummy),
	(r"\\subsection{(.*?)}", (lambda : r"===\1==="), dummy),
	(r"\\section{(.*?)}", (lambda : r"==\1=="), dummy),
	(r"\\_", (lambda :"_"), dummy),
#	(r"\\title{(.*)}", (lambda :r"= \1 ="),dummy),
#        (r"\\author{(.*)}", (lambda :r"\1"),dummy),
        (r"\\date{(.*)}", (lambda :r"\1"),dummy),
	(r"\\tableofcontents",None, dummy),
	(r"\\null",None, dummy),
	(r"\\newpage",None, dummy),
	(r"\\thispagestyle{.*?}", None, dummy),
	(r"\\maketitle", None, dummy),
	(r"\n$", decide_el, dummy),
#	(r"[^\\]?\{", None, dummy),
#	(r"[^\\]?\}", None, dummy),
        (r"\\\$",(lambda : r"&#36;"),dummy),
#	(r"\$(.*?)\$",(lambda :r"<math>\1</math>"),dummy),
	(r"\$",decide_math,toggle_math),
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
        (r"\\LaTeX",(lambda: "L<sup>A</sup>TEX"), dummy),
        (r"\\\\", None, dummy),
        (r"\\ ",(lambda: " "), dummy),
        # unknown command
        (r"\\\w+{[^\}]+}",None, dummy),
]

# precompile regular expressions
reg = map(lambda x: (re.compile(x[0]),x[1],x[2]),tr_list2)

def convert(in_stream,out_stream):
    "Convert LaTeX to MediaWiki"
    for i in in_stream.readlines():
	mystr = i

	for r in reg:
            if r[0].search(mystr):
                r[2]()
            if r[1]:
                mysub = r[1]()
            else:
                mysub = ""                
            mystr = r[0].sub(mysub, mystr)

	if bdoc:
            out_stream.write(mystr)

if __name__ == "__main__":
    convert(sys.stdin,sys.stdout)
