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

pclip=par.float("pclip",100) # clip percent

# style
plotfat=par.int("plotfat", 1) # plot line width
plotcol=par.string("plotcol", "0,1,2,3,4,5,6,7") 
# plot color
dash=par.string("dash", "0") 
# dash styles \n 0	solid line \n 1	dash line \n 2	dotted line \n 3	dash dot
symbol=par.string("symbol", None) # mark symbols
legends=par.string("legend", ",") # legends
plotcol=plotcol.split(',')
dash=dash.split(',')
legends=legends.split(',')
if symbol != None:
	symbol=symbol.split(',')


x1=arange(n1)*d1+o1
x=zeros((n2,n1),'f')
input.read(x)

xx=(x.max()-x.min())*(100-pclip)/100
min2=x.min()+xx
max2=x.max()-xx

min1=par.float("min1",x1[0])
max1=par.float("max1",x1[n1-1])
min2=par.float("min2",min2)
max2=par.float("max2",max2)

font = {'family' : 'serif'}

rc('font', **font)

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
label1=par.string("label1")
if label1 == None:
	label1=input.string("label1")
if label1 != None:
	label=label1
	unit1=input.string("unit1")
	if unit1 == None:
		unit1=par.string("unit1",None)
	if unit1 != None:
		label =label+': ('+ unit1+')'
	xlabel(label)

label2=par.string("label2")
if label2 == None:
	label2=input.string("label2")
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

xlim(min1,max1)
ylim(min2,max2)


savefig(sys.stdout, format='pdf', bbox_inches='tight', transparents=True)

