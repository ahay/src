#!/usr/bin/env python
##   Copyright (C) 2008 University of Texas at Austin
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
import math
import vplot

def draw(vp,NARRAY=41):
    plotfat=1
    wantframe=False
    option=True
    animate=False
    xll=1.705
    yll=1.37
    xur=11.945
    yur=8.87
    h=1.0
    y=0.0
    time=1.6
    v=2.0
    theta=30

    a = v * time /2
    b = math.sqrt(a*a-h*h)
    xs = y - h
    xg = y + h
    surface = 0

    slop = 0.05*a
    xmin = -2.
    zmin = -1.
    xmax =  2.
    zmax = 2.
    xcenter=0.5*(xmin+xmax)
    zcenter=0.5*(zmin+zmax)

    # set vplot scales 
    vp.orig(xll,yur)
    vp.uorig(xmin,zmin)
    xscale = (xur-xll)/(xmax-xmin)
    zscale = -(yur-yll)/(zmax-zmin)
    vp.scale(xscale,zscale)

    # 					draw and annotate surface 
    vp.fat(2*plotfat)

    vp.color(7)   
    vp.umove(xmin,surface)
    vp.udraw(xmax,surface)

    vp.color(5)
    vp.fat(plotfat)
    vp.tjust(4,6)
    vp.utext(xs-0.5*slop,surface-slop,8,0,"S")
    vp.utext(xg,surface-slop,8,0,"G")
    vp.utext(y,surface-slop,8,0,"M")

    dth = 180.0/(NARRAY-1)

    vp.color(4)
    # LOOP on theta (ith) to generate ellipse
    xval = []
    zval = []
    for ith in xrange(NARRAY):
        thtmp = 90. - ith * dth
        thtmp = thtmp*math.pi/180.
        tanth = math.tan(thtmp)
        sinth = math.sin(thtmp)
	costh = math.cos(thtmp)
        
	d = math.hypot(a*sinth,b*costh)

	xr  = y - d*sinth - h*h*costh*costh*sinth/d
	zr  = surface + d*costh - h*h*costh*sinth*sinth/d
        
	xval.append(xr)
	zval.append(zr)
	
    #					Draw ellipse 
    vp.color(5)
    vp.upline(xval,zval,NARRAY)

    if option:
        #	  Loop on theta to draw SRG raypaths 

	iflag = 0
	emph  = 3

        for ith in xrange(0,NARRAY,5):
	    iflag += 1
	    thtmp = 90. - ith * dth
	    thtmp = thtmp*math.pi/180
	    tanth = math.tan(thtmp)
	    sinth = math.sin(thtmp)
	    costh = math.cos(thtmp)

	    d = math.hypot(a*sinth,b*costh)

	    xr  = y - d*sinth - h*h*costh*costh*sinth/d
	    zr  = surface + d*costh - h*h*costh*sinth*sinth/d

            #					Shot to reflection point 
	    vp.color(4)
            
	    x = xs
	    z = surface
	    x1 = xr
	    z1 = zr
	    if iflag == emph:
                vp.color(6)
	    if ith != 1 and ith != NARRAY:
                vp.uarrow(x,z,x1,z1,0.02*zmax)
	    if iflag == emph:
                vp.color(4)
            #					Annotate beta 
	    xend = (xs + xr)/2.
	    zend = (surface + zr)/2.

	    if iflag == emph:
		vp.color(6)
		vp.utext(xend-slop,  zend+0.8*slop,8,0,"\\F9 b")
		vp.color(4)
                

            #					reflection point to geophone 
	    x = xr
	    z = zr
	    vp.umove(x,z)
	    x1 = xg
	    z1 = surface
	    if iflag == emph:
                vp.color(6)

	    if ith != 1 and ith != NARRAY:
                vp.uarrow(x,z,x1,z1,0.02*zmax)

	    if iflag == emph:
                vp.color(6)
            #					Annotate alpha 
	    xend = (xg + xr)/2.0
	    zend = (surface + zr)/2.0

	    if iflag == emph:
		vp.color(6)
		vp.utext(xend+slop,  zend+0.5*slop,8,0,"\\F9 a")
		vp.color(4)
	    
            if animate:
		vp.erase()	      # clear screen, then redraw basic stuff 
                # 			redraw and annotate surface 
		vp.fat(2*plotfat)
		vp.color(7)
		vp.umove(xmin,surface)		
		vp.udraw(xmax,surface)

		vp.color(5)
		vp.fat(plotfat)
		vp.tjust(4,6)
		vp.utext(xs-0.5*slop,surface-slop,8,0,"S")
		vp.utext(xg,surface-slop,8,0,"G")
		vp.utext(y,surface-slop,8,0,"M")

                #					redraw ellipse 
		vp.color(5)
		vp.upline(xval,zval,NARRAY)
    else:
        #	                  Begin work for drawing dipping layer 
        theta *= math.pi/180
        tanth = math.tan(theta)
        sinth = math.sin(theta)
	costh = math.cos(theta)

	d = math.hypot(a*sinth,b*costh)

	xr  = y - d*sinth - h*h*costh*costh*sinth/d
	zr  = surface + d*costh - h*h*costh*sinth*sinth/d

	xz  = y - h*h/d*sinth
	xsp = y - h   - 2*(d-h*sinth)*sinth
        xrp = y - d*sinth
        zrp = surface + d*costh
	xstart = xsp
	zstart = surface + d/costh + tanth * (xstart-y)
	xend = 0.5*(xg + y)
	zend = surface + d/costh + tanth * (xend-y)


        # 					Draw dipping reflector 
	vp.color(2)
	vp.umove(xstart,zstart)
	vp.udraw(xend,zend)
	vp.dash(0.,0.,0.,0.)
        #					indicate angle theta 
	vp.color(7)
	vp.dash(.1,.1,.1,.1)
	vp.umove(xend,zend)
	vp.udraw(xend-0.5*h,zend)
	vp.dash(.0,.0,.0,.0)

        # 					finite-offset raypaths 
	vp.color(4)

        #					Shot to reflection point 
 	vp.uarrow(xs,surface,xr,zr,0.02*zmax)
	
        #					reflection point to geophone 
	vp.uarrow(xr,zr,xg,surface,0.02*zmax)
	
        # 					text 


        vp.color(5)
        vp.fat(plotfat)
        vp.tjust(4,6)
        vp.utext(xr,zr+2.0*slop,8,0,"R")
 	vp.utext(xend-5*slop,  zend-1.*slop,8,0,"\\F9 q")
    

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
    draw(vp)


