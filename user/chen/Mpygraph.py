#!/usr/bin/env python
'graph by python'

##   Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
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

import sys
try:
	from pylab import *
	from string import *
	import rsf.api as rsf
except Exception, e:
	print 'ERROR: pylab needed'
	sys.exit(1)

colortable=('black','blue','red','magenta','green','cyan','yellow','white')
styletable=('-', '--', ':' , '-.')

par=rsf.Par()
input=rsf.Input()
n1=input.int("n1")
o1=input.float("o1")
d1=input.float("d1")
n2=input.int("n2")
o2=input.float("o2")
d2=input.float("d2")

pclip=par.float("pclip",100) # clip percent

# style
plotfat=par.int("plotfat", 1) # plot line width
markint=par.int("markint", 10) # mark symbol interval
plotcol=par.string("plotcol", "0,1,2,3,4,5,6,7") 
# plot color
dash=par.string("dash", "0") 
# dash styles \n 0	solid line \n 1	dash line \n 2	dotted line \n 3	dash dot
symbol=par.string("symbol", None) # mark symbols
legends=par.string("legend", None) # legends
xpos=par.string("x", None) # xpos
plotcol=plotcol.split(',')
dash=dash.split(',')
if legends != None:
	legends=legends.split(':')
if symbol != None:
	symbol=symbol.split(',')

usetex=par.bool("usetex", False) # use tex symbol


x=zeros((n2,n1),'f')
input.read(x)
if xpos != None:
	xfile = rsf.Input(xpos)
	x1=zeros(n1,'f')
	xfile.read(x1)
else:
	x1=arange(n1)*d1+o1


xx=(x.max()-x.min())*(100-pclip)/100
min2=x.min()+xx
max2=x.max()-xx

min1=par.float("min1",x1[0])
max1=par.float("max1",x1[n1-1])
min2=par.float("min2",min2)
max2=par.float("max2",max2)

font = {'family' : 'serif'}

rc('text', usetex=usetex)
rc('font', **font)

for i2 in range (n2):
	if legends != None:
		if i2 < len(legends):
			label0=legends[i2]
		else:
			label0=None
	else:
		label0=None
	id=int(plotcol[i2%len(plotcol)])%len(colortable)
	color0=colortable[id]
	id=int(dash[i2%len(dash)])%len(styletable)
	ls0=styletable[id]
	if symbol != None and i2 < len(symbol):
		marker0=symbol[i2]
	else:
		marker0='None'
	plot(x1,x[i2,:],
	lw=plotfat, ls=ls0,
	marker=marker0, markevery=n1/markint,
	color=color0,
	label=label0)


# x label
label1=par.string("label1")
if label1 == None:
	label1=input.string("label1")
if label1 != None:
	label=label1
	unit1=input.string("unit1")
	if unit1 == None:
		unit1=par.string("unit1",None)
	if unit1 != None:
		label =label+' ('+ unit1+')'
	ylabel(label)

label2=par.string("label2")
if label2 == None:
	label2=input.string("label2")
if label2 != None:
	label=label2
	unit2=input.string("unit2")
	if unit2 == None:
		unit2=par.string("unit2",None)
	if unit2 != None:
		label =label+' ('+ unit2+')'
	xlabel(label)


# legends:
wantlegend=par.bool("wantlegend", True)
wherelegend=par.int("wherelegend", 1)
if wantlegend:
	legend(loc=wherelegend)

xlim(min1,max1)
ylim(min2,max2)

format=par.string("format","pdf")

savefig(sys.stdout, format=format, bbox_inches='tight', transparents=True)

