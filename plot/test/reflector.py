#!/usr/bin/env python
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

import sys, math

import vplot

def draw(vp,what,wantframe=False):
    narcs=3
    xll=1.705
    yll=1.37
    xur=11.945
    yur=8.87
    plotfat=1
    NANGLES=11
    NARCS2=6
    d=1.4
   
    surface=0
    y=0.5
    theta=30
    
    if what=='expt':
        arcln=45
    else:
        arcln=30 
    
    degrad = math.pi/180
    theta *= degrad;
    tanth = math.tan(theta)
    sinth = math.sin(theta)
    costh = math.cos(theta)

    dth = arcln*degrad/(NANGLES-1)
    arcmin = -arcln / 2. * degrad
    drad = d / (narcs + 1)

    xs = y					# Shotpoint 
    zs = surface
    xr = xs - d*sinth 			        # Reflection point 
    zr = zs + d*costh 
    xc = xr					# Migrated location 
    zc = surface
    xsp = y  - 2*d*sinth			# Image shotpoint 
    zsp = surface + 2*d*costh

    xs1= y + d/2.				# Second shotpoint 
    zs1= surface
    d1 = d + (xs1-xs)*sinth			# Second Reflection point 
    xr1= xs1 - d1*sinth			
    zr1= zs1 + d1*costh 

    xs2= y - d/2.				# Third shotpoint 
    zs2= surface
    d2 = d + (xs2-xs)*sinth			# Third Reflection point 
    xr2= xs2 - d2*sinth			
    zr2= zs2 + d2*costh 

    slop = 0.05*d
    xmin = -2.
    zmin = 0.
    xmax =  2.
    zmax = 3.
    xcenter=0.5*(xmin+xmax)
    zcenter=0.5*(zmin+zmax)

    #						dipping-reflector parameters 
    xstart = xc - d
    if  xstart < xmin:
        xstart = xmin
    xend = xc + d
    if  xend > xmax: 
        xend = xmax
    zstart = surface + d/costh + tanth * (xstart-y)
    zend   = surface + d/costh + tanth * (xend-y)

#						Pre-compute arcs of circles 

#						Downgoing wavefronts 
    radius = 0.
    xval = []
    zval = []
    for iarc in range(narcs):
	radius += drad
        xvali = []
        zvali = []
	for iang in range(NANGLES):
	    phi = arcmin + iang * dth
	    xvali.append(xs - radius * math.sin(theta+phi))
	    zvali.append(zs + radius * math.cos(theta+phi))
	xval.append(xvali)
        zval.append(zvali)

#						Upcoming wavefronts 
    radius = d + drad/2
    for iarc in range(NARCS2-narcs):
	radius += drad
        xvali = []
        zvali = []
        for iang in range(NANGLES):
	    phi = arcmin + iang * dth
	    xvali.append(xsp + radius * math.sin(theta+phi))
	    zvali.append(zsp - radius * math.cos(theta+phi))
        xval.append(xvali)
        zval.append(zvali)  

    vp.orig(xll,yur)
    vp.uorig(xmin,zmin)
    xscale = (xur-xll)/(xmax-xmin)
    zscale = -(yur-yll)/(zmax-zmin)
    vp.scale(xscale,zscale)
#						Text centering
    vp.tjust(4,1)

# 					draw surface 
    vp.fat(2*plotfat)

    vp.color(7)
    vp.umove(xmin,surface)
    vp.udraw(xmax,surface)

#					draw dipping reflector 
    vp.color(2)   
    vp.umove(xstart,zstart)
    vp.udraw(xend,zend)
    vp.dash(0.,0.,0.,0.)

#					indicate angle theta 
    vp.color(5)
    vp.utext(xend-5*slop,zend-0.5*slop,8,0,"\\F9 q")

    vp.color(7)
    vp.dash(.1,.1,.1,.1)
    vp.umove(xend,zend)
    vp.udraw(xend-0.5*d,zend)
    vp.dash(.0,.0,.0,.0)

# 					main zero-offset raypath
    vp.color(4)
#					Shot to reflection point
    vp.uarrow(xs,surface,xr,zr,0.02*zmax)	
#					reflection point to geophone
    vp.uarrow(xr,zr,xs,surface,0.02*zmax)
	
    if what=='expt':
        vp.color(5)
        vp.fat(plotfat)
        vp.utext(xs, zs -slop,8,0,"S\\s60 \\_2")
        vp.utext(xs1,zs1-slop,8,0,"S\\s60 \\_3")
        vp.utext(xs2,zs2-slop,8,0,"S\\s60 \\_1")
        vp.utext(xr,zr+3.0*slop,8,0,"R\\s60 \\_2")
        vp.utext(xr1,zr1+3.0*slop,8,0,"R\\s60 \\_3")
        vp.utext(xr2,zr2+3.0*slop,8,0,"R\\s60 \\_1")
    #					  second shotpoint 
        vp.color(4)
        vp.uarrow(xs1,zs1,xr1,zr1,0.02*zmax)
        vp.uarrow(xr1,zr1,xs1,zs1,0.02*zmax)

    #					  third shotpoint 
        vp.color(4)
        vp.uarrow(xs2,zs2,xr2,zr2,0.02*zmax)
        vp.uarrow(xr2,zr2,xs2,zs2,0.02*zmax)

    #					  Draw Wavefronts 
        vp.color(6)
        for iarc in range(NARCS2):
            vp.upline(xval[iarc],zval[iarc],NANGLES)
    else:			      
        # reflector kinematics 
        vp.color(5)
        vp.fat(plotfat)
        vp.utext(xs,surface-slop,8,0,"S")
        vp.utext(xc,surface-slop,8,0,"C")
        vp.utext(xr,zr+3.0*slop,8,0,"R")
        x=0.5*(xr+y)
        z=0.5*(zr+surface)
        vp.utext(x-slop,z,8,0,"d")
        vp.dash(.1,.1,.1,.1)	# dashed vertical line 
        vp.umove(xc,zc)
        vp.udraw(xr,zr)
        x=xc 
        z=0.5*(zc+zr)
        vp.utext(x-slop,z,8,0,"z")

    if wantframe:
	vp.dash(0.,0.,0.,0.)
	vp.color(7)
	vp.fat(1)
	vp.umove(xmin,zmin)
	vp.udraw(xmin,zmax)
	vp.udraw(xmax,zmax)
	vp.udraw(xmax,zmin)
	vp.udraw(xmin,zmin)

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp,sys.argv[1])
