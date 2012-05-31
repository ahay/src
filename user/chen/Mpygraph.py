#!/usr/bin/env python
'graph by python'

import sys
try:
	from pylab import *
	from string import *
	import rsf.api as rsf
except Exception, e:
	print 'ERROR: numpy needed'
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

pclip=par.float("pclip",100)

# style
plotfat=par.int("plotfat", 1)
plotcol=par.string("plotcol", None)
dash=par.string("dash", None)
symbol=par.string("symbol", None)
legends=par.string("legend", None)
if plotcol != None:
	plotcol=plotcol.split(',')
else:
	plotcol=arange(n2)
if dash != None:
	dash=dash.split(',')
else:
	dash=zeros(n2,'int')
if legends != None:
	legends=legends.split(',')
else:
	legends=arange(n2)
if symbol != None:
	symbol=symbol.split(',')

#plotcol=plotcol%len(colortable)
#dash=dash%len(styletable)

#sys.stderr.write(leg2.tostring())

x1=arange(n1)*d1+o1
x=zeros((n2,n1),'f')
input.read(x)

xx=(x.max()-x.min())*(100-pclip)/100
min=x.min()+xx
max=x.max()-xx


for i2 in range (n2):
	if i2 < len(legends):
		label0=legends[i2]
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
	marker=marker0, markevery=n1/10,
	color=color0,
	label=label0)


# x label
label1=input.string("label1")
if label1 == None:
	label1=par.string("label1")
if label1 != None:
	label=label1
	unit1=input.string("unit1")
	if unit1 == None:
		unit1=par.string("unit1",None)
	if unit1 != None:
		label =label+': ('+ unit1+')'
	xlabel(label)

label2=input.string("label2")
if label2 == None:
	label2=par.string("label2")
if label2 != None:
	label=label2
	unit2=input.string("unit2")
	if unit2 == None:
		unit2=par.string("unit2",None)
	if unit2 != None:
		label =label+': ('+ unit2+')'
	ylabel(label)


# legends:
wantlegend=par.bool("wantlegend", True)
wherelegend=par.int("wherelegend", 1)
if wantlegend:
	legend(loc=wherelegend)

if pclip < 100 and pclip >0:
	ylim(min,max)


savefig(sys.stdout, format='pdf')

