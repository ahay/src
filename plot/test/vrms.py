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
import sys, math
import vplot

# quantities related to vplot scales
zbottom = 0.
dz=(1.0,2.0,1.5)
for h in dz:
    zbottom += h
    
theta3 = 20*math.pi/180
sn=[0,0,math.sin(theta3)]
cs=[0,0,0]
vv=[1.5,2.0,2.5]

h=0.
t=0.
for i in (2,1,0):
    if i != 2:
        sn[i] = sn[i+1]*vv[i]/vv[i+1] # snell's law
    cs[i] = math.sqrt(1. - sn[i]*sn[i])
    h += dz[i]*sn[i]/cs[i]
    t += 2*dz[i]/(cs[i]*vv[i])


xmin = -h - 0.1*zbottom
if  xmin > -2:
    xmin = -2

zmin = -0.1 * zbottom
xmax = -xmin
zmax = 1.1*zbottom
if  zmax < 3:
    zmax = 3

if __name__ == "__main__":
    vp = vplot.Vplot()

    xll=1.705
    yll=1.37
    xur=11.945
    yur=8.87

    vp.orig(xll,yur)
    vp.uorig(xmin,zmin)
    
    xscale = (xur-xll)/(xmax-xmin);
    zscale = -(yur-yll)/(zmax-zmin);
    vp.scale(xscale,zscale);

    # draw flat reflectors 
    vp.color(6)
    vp.fat(1)

    z=0.
    for d in dz:
	z += d
	vp.umove(xmin,z)
	vp.udraw(xmax,z)
    

    # draw vertical zero-offset ray
    vp.color(7)
    vp.umove(0.,0.)
    vp.udraw(0.,zbottom)

    # draw ray paths from bottom reflector 
    vp.color(2)
    rad = 0.02*zbottom

    for sign in (1,-1):
        x=0 
        z=zbottom
        for i in (2,1,0):
            if i != (1+sign):
                vp.umove(x,z)
            else:
                x1 = x 
                z1 = z
            
            x -= sign*dz[i]*sn[i]/cs[i]
            z -= dz[i]

            if i != (1+sign):
                vp.udraw(x,z)
            else:
                if sign==1:
                    vp.uarrow(x,z,x1,z1,rad)
                else:
                    vp.uarrow(x1,z1,x,z,rad)

    # draw surface
    vp.color(7)
    vp.fat(2)
    vp.umove(xmin,0)
    vp.udraw(xmax,0)

    
    # draw text
    vp.color(5)
    vp.fat(2)
    vp.tjust(2,3)

    x = 0.
    z = zbottom * 1.05
    vp.utext(x,z,8,0,"R")

    x = -h
    z = 0.- zbottom *.05
    vp.utext(x,z,8,0,"S")

    x = h
    z = 0. - zbottom *.05
    vp.utext(x,z,8,0,"G")

    x = -h
    z = 0. + dz[0]/2.
    vp.utext(x,z,8,0,"v\\s60 \\_1")

    z += dz[0]/2. + dz[1]/2.
    vp.utext(x,z,8,0,"v\\s60 \\_2")

    z += dz[1]/2. + dz[2]/2.
    vp.utext(x,z,8,0,"v\\s60 \\_3")

    x = 0
    z = 0. - zbottom *.05
    vp.utext(x,z,8,0,"2h")
    x = 0. + 0.05*zbottom
    vp.umove(x,z)
    x1 = h - 0.05 * zbottom
    vp.uarrow(x,z,x1,z,-0.02*zbottom)
    x = 0. - 0.05*zbottom
    vp.umove(x,z)
    x1 = -h + 0.05 * zbottom
    vp.uarrow(x,z,x1,z,-0.02*zbottom)

    # Highlight the trigonometry
    vp.color(7)
    vp.fat(1)

    # bottom of triangle
    z=zbottom-dz[2]
    x=dz[2]*sn[2]/cs[2]
    vp.umove(x,z)
    x1 = x + dz[1]*sn[1]/cs[1]
    z1 = z
    vp.udraw(x1,z1)
    xtxt = (x+x1)/2.
    ztxt = z + 0.05 * zbottom
    vp.utext(xtxt,ztxt,8,0,"\\F9 D\\F3 x\\_\\s60 i")

    # vertical side of triangle
    x = x1
    z = z1 - dz[1]
    vp.umove(x1,z1)
    vp.udraw(x,z)
    xtxt = x + 0.05 * zbottom
    ztxt = (z+z1)/2.
    vp.utext(xtxt,ztxt,8,0,"\\F9 D\\F3 z\\_\\s60 i")
    xtxt = x - sn[1]/cs[1]*dz[1]/8.
    ztxt = z + dz[1]/4.
    vp.utext(xtxt,ztxt,8,0,"\\F9 q\\F3 \\s60 \\_i")

    # vertical side in tau space 
    xtxt = 0. - 0.05*zbottom
    ztxt = z + dz[1]/2.
    vp.utext(xtxt,ztxt,8,90,
	     "\\F9 Dt\\F3 \\s60 \\_i\\^\\s100 v\\s60 \\_i\\^\\s100 /2")

    # hypotenuse of triangle
    z1=zbottom-dz[2]
    x1=dz[2]*sn[2]/cs[2]
    vp.umove(x,z)
    vp.udraw(x1,z1)
    xtxt = (x+x1)/2. - 0.05 *zbottom
    ztxt = (z+z1)/2.
    vp.utext(xtxt,ztxt,8,60,
	     "\\F9 D\\F3 t\\s60 \\_i\\^\\s100 v\\s60 \\_i\\^\\s100 /2")

    wantframe=False
    if wantframe:
	vp.umove(xmin,zmin)
	vp.udraw(xmin,zmax)
	vp.udraw(xmax,zmax)
	vp.udraw(xmax,zmin)
	vp.udraw(xmin,zmin)
