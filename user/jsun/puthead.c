/*************************************************************************
* *  put head for rsf datasets
* *************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
             
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
                   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

void puthead1(sf_file Fo, int n1, int n2, float d1, float o1, float d2, float o2)
/*<  put head for (kx,kz) domain float-type data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","kz");
        sf_putstring(Fo,"label2","kx");
        sf_putstring(Fo,"unit1","2*pi/m");
        sf_putstring(Fo,"unit2","2*pi/m");
}

void puthead2(sf_file Fo, int n1, int n2, float d1, float d2, float o1, float o2)
/*<  put head for (x,z) domain float-type 2D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"unit1","m");
        sf_putstring(Fo,"unit2","m");
}

void puthead2x(sf_file Fo, int n1, int n2, float d1, float d2, float o1, float o2)
/*<  put head for (x,t) domain float-type 2D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","t");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"unit1","second");
        sf_putstring(Fo,"unit2","m");
}

void puthead2kx(sf_file Fo, int n1, int n2, float d1, float d2, float o1, float o2)
/*<  put head for (kx,kz) domain float-type 2D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","kz");
        sf_putstring(Fo,"label2","kx");
        sf_putstring(Fo,"unit1","2*pi/m");
        sf_putstring(Fo,"unit2","2*pi/m");
}

void puthead3(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (z,x,t) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putfloat(Fo,"o3",o3);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"label3","t");
        sf_putstring(Fo,"unit1","m");
        sf_putstring(Fo,"unit2","m");
        sf_putstring(Fo,"unit3","second");
}

void puthead3x(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (t,x,s) shot gather cubes  >*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putfloat(Fo,"o3",o3);
        sf_putstring(Fo,"label1","t");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"label3","shot");
        sf_putstring(Fo,"unit1","second");
        sf_putstring(Fo,"unit2","m");
        sf_putstring(Fo,"unit3","m");
}

void puthead3kx(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (ky,kx,kz) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putfloat(Fo,"o3",o3);
        sf_putstring(Fo,"label1","kz");
        sf_putstring(Fo,"label2","kx");
        sf_putstring(Fo,"label3","ky");
        sf_putstring(Fo,"unit1","2*pi/m");
        sf_putstring(Fo,"unit2","2*pi/m");
        sf_putstring(Fo,"unit3","2*pi/m");
}

void puthead2dcommonshot(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for 2D (x,t) domain common-shot data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putfloat(Fo,"o3",o3);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"label3","y");
        sf_putstring(Fo,"unit1","km");
        sf_putstring(Fo,"unit2","km");
        sf_putstring(Fo,"unit3","km");
}
