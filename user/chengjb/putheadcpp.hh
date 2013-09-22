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


#include <rsf.hh>

static void puthead1(oRSF Fo, int n1, int n2, float d1, float o1, float d2, float o2)
/*<  put head for (kx,kz) domain float-type data sets>*/
{
        /* Read/Write axes */
        Fo.put("n1",n1);
        Fo.put("n2",n2);
        Fo.put("d1",d1);
        Fo.put("o1",o1);
        Fo.put("d2",d2);
        Fo.put("o2",o2);
        Fo.put("label1","kz");
        Fo.put("label2","kx");
        Fo.put("unit1","2*pi/m");
        Fo.put("unit2","2*pi/m");
}

static void puthead2(oRSF Fo, int n1, int n2, float d1, float o1, float d2, float o2)
/*<  put head for (x,z) domain float-type data sets>*/
{
        /* Read/Write axes */
        Fo.put("n1",n1);
        Fo.put("n2",n2);
        Fo.put("d1",d1);
        Fo.put("o1",o1);
        Fo.put("d2",d2);
        Fo.put("o2",o2);
        Fo.put("label1","z");
        Fo.put("label2","x");
        Fo.put("unit1","km");
        Fo.put("unit2","km");
}

static void puthead3(oRSF Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (x,z,t) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        Fo.put("n1",n1);
        Fo.put("n2",n2);
        Fo.put("n3",n3);
        Fo.put("d1",d1);
        Fo.put("d2",d2);
        Fo.put("d3",d3);
        Fo.put("o1",o1);
        Fo.put("o2",o2);
        Fo.put("o3",o3);
        Fo.put("label1","z");
        Fo.put("label2","x");
        Fo.put("label3","t");
        Fo.put("unit1","km");
        Fo.put("unit2","km");
        Fo.put("unit3","second");
}

static void puthead3x(oRSF Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (y,x,z) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        Fo.put("n1",n1);
        Fo.put("n2",n2);
        Fo.put("n3",n3);
        Fo.put("d1",d1);
        Fo.put("d2",d2);
        Fo.put("d3",d3);
        Fo.put("o1",o1);
        Fo.put("o2",o2);
        Fo.put("o3",o3);
        Fo.put("label1","z");
        Fo.put("label2","x");
        Fo.put("label3","y");
        Fo.put("unit1","km");
        Fo.put("unit2","km");
        Fo.put("unit3","km");
}

static void puthead3kx(oRSF Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2, float o3)
/*<  put head for (ky,kx,kz) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        Fo.put("n1",n1);
        Fo.put("n2",n2);
        Fo.put("n3",n3);
        Fo.put("d1",d1);
        Fo.put("d2",d2);
        Fo.put("d3",d3);
        Fo.put("o1",o1);
        Fo.put("o2",o2);
        Fo.put("o3",o3);
        Fo.put("label1","kz");
        Fo.put("label2","kx");
        Fo.put("label3","ky");
        Fo.put("unit1","2*pi/m");
        Fo.put("unit2","2*pi/m");
        Fo.put("unit3","2*pi/m");
}
