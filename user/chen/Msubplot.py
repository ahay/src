#!/usr/bin/env python
'subplot by python'

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

def CustomLocator(ot, dt, min, max):
	ticks=[]
	tic=ot
	if dt>0:
		while (tic < max):
			ticks.append(tic)
			tic=tic+dt
	if dt<0:
		while (tic > min):
			ticks.append(tic)
			tic=tic+dt
	return FixedLocator(ticks)

def TicLocator(ntic, otic, dtic, min, max):
	if ntic == 0:
		tic=AutoLocator()
	elif ntic==1:
		tic=CustomLocator(otic, dtic, min, max)
	elif ntic < 0 :
		tic=NullLocator()
	else :
		tic=MultipleLocator(ntic)
	return tic

	

def TicFormatter(id):
	if id==1:
		fmt=FormatStrFormatter('%1.1f')
	if id==2:
		fmt=FormatStrFormatter('%1.2f')
	if id <0:
		fmt=ScalarFormatter(useOffset=False,useMathText=True)
		fmt.set_powerlimits((id,id+1))
		fmt.set_scientific(True)
	else:
		fmt=ScalarFormatter()
	return fmt

def getnstr(par,key,default, nx=1):
	xstr=par.string(key)
	if(xstr == None):
		xx=[]
		xx.append(default)
	else:
		xx=xstr.split(',')
		if xx == None:
			xx.append(default)
	while (len(xx)<nx):
		xx.append(default)
	return xx

def getnint(par,key,default, nx=1, min=0, max=10000):
	xstr=par.string(key)
	xx=[]
	if(xstr == None):
		xx.append(default)
	else:
		xlst=xstr.split(',')
		if(xlst == None):
			xx.append(default)
		for ii in xlst:
			x0=int(ii)
			if x0<min and x0>max:
				xx.append(default)
			else :
				xx.append(x0)
	while (len(xx)<nx):
		xx.append(default)
	return xx


def getnfloat(par,key,default, nx=1, min=0, max=10000):
	xstr=par.string(key)
	xx=[]
	if(xstr == None):
		xx.append(default)
	else:
		xlst=xstr.split(',')
		if(xlst == None):
			xx.append(default)
		for ii in xlst:
			x0=float(ii)
			if x0<min and x0>max:
				xx.append(default)
			else :
				xx.append(x0)
	while (len(xx)<nx):
		xx.append(default)
	return xx

colortable=('black','blue','red','magenta','green','cyan','yellow','white')
styletable=('-', '--', ':' , '-.')
font = {'family' : 'serif'}
rc('font', **font)

par=rsf.Par()
input=rsf.Input()
n1=input.int("n1")
o1=input.float("o1")
d1=input.float("d1")
n2=input.int("n2")
o2=input.float("o2")
d2=input.float("d2")


np1=par.int("np1",1) # subplot(np2,np1)
np2=par.int("np2",n2) # subplot(np2,np1)

np=np1*np2

x1=arange(n1)*d1+o1
x=zeros((n2,n1),'f')
input.read(x)

label1=input.string("label1")
label2=input.string("label2")
unit1=input.string("unit1")
unit2=input.string("unit2")
if label1 == None:
	label1=""
if label2 == None:
	label2=""
if unit1 == None:
	unit1=""
if unit2 == None:
	unit2=""
label1=getnstr(par,"label1",label1, nx=np)
label2=getnstr(par,"label2",label2, nx=np)
unit1=getnstr(par,"unit1",unit1, nx=np)
unit2=getnstr(par,"unit2",unit2, nx=np)


scale=getnint(par,"scale", 0, nx=np, min=-5, max=5)
fmt1=getnint(par,"fmt1", 0, nx=np, min=-5, max=5)
fmt2=getnint(par,"fmt2", 0, nx=np, min=-5, max=5)
n1tic=getnint(par,"n1tic", 0, nx=np, min=-10000)
n2tic=getnint(par,"n2tic", 0, nx=np, min=-10000)
o1tic=getnfloat(par, "o1tic", 0.0, nx=np, min=-10000.0)
d1tic=getnfloat(par, "d1tic", 1.0, nx=np, min=-10000.0)
o2tic=getnfloat(par, "o2tic", 0.0, nx=np, min=-10000.0)
d2tic=getnfloat(par, "d2tic", 1.0, nx=np, min=-10000.0)

for ip in range(np):
	dat=getnint(par, 'p%d'%ip, ip, nx=-1, max=n2)
	xx=x[dat,:]*10**scale[ip]
	nd=len(dat)
	fat=getnint(par, 'p%dfat'%ip, 1, nx=nd, max=5) # plot line width
	col=getnint(par, 'p%dcol'%ip, 0, nx=nd, max=7) # plot line width
	dash=getnint(par, 'p%ddash'%ip, 0, nx=nd, max=3) # plot line width
	lgd=getnstr(par, 'p%dlegend'%ip, "", nx=nd) # legends
	syb=getnstr(par, 'p%dsymbol'%ip, "", nx=nd) # mark symbols
	pclip=par.float("p%dclip"%ip, 0.0)
	min1=par.float("p%dmin1"%ip, x1.min())
	max1=par.float("p%dmax1"%ip, x1.max())
	ll=(xx.max()-xx.min())*pclip
	min2=par.float("p%dmin2"%ip, xx.min()+ll)
	max2=par.float("p%dmax2"%ip, xx.max()-ll)

	ax=subplot(np2,np1,ip+1)
	ax.xaxis.set_major_formatter(
		TicFormatter(fmt1[ip]))
	ax.yaxis.set_major_formatter(
		TicFormatter(fmt2[ip]))
	ax.xaxis.set_major_locator(
		TicLocator(n1tic[ip], o1tic[ip], d1tic[ip], min1, max1))
	ax.yaxis.set_major_locator(
		TicLocator(n2tic[ip], o2tic[ip], d2tic[ip], min2, max2))
	for id in range(nd):
		ax.plot(x1,xx[id,:],
			color=colortable[col[id]],
			lw=fat[id], 
			marker=syb[id], markevery=n1/10,
			ls=styletable[dash[id]])
	if label1[ip] != "":
		label=label1[ip]
		if unit1[ip] != "":
			label='%s (%s)'%(label,unit1[ip])
		xlabel(label)
	if label2[ip] != "":
		label=label2[ip]
		if unit2[ip] != "":
			label='%s (%s)'%(label,unit2[ip])
		ylabel(label)
	xlim(min1,max1)
	ylim(min2,max2)
	if(scale[ip]!=0):
		locs,labs=yticks()
		text(min1+0.01*(max1-min1), min2+0.8*(max2-min2),
			r'$\times 10^{%d}$'%(-scale[ip]))


savefig(sys.stdout, format='pdf', bbox_inches='tight', transparents=True)

