/* 2-D two-components wavefield modeling using pseudo-pure mode P-wave equation in VTI media.
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng and Wei Kang
     
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

/* prepared head files by myself */
#include "_cjb.h"

/* calculate phase, group velocity and angle  */
#include "zero.h"

int main(int argc, char* argv[])
{

        int    nx, nz, na, nb, i, j;
		float  fx, fz, dx, dz, da, db;

        sf_init(argc,argv);

        sf_file Fo1, Fo2, Fo3, Fo4, Fo5, Fo6;
        Fo1 = sf_input("in"); // wavefront1
        Fo2 = sf_input("wave2"); // wavefront2
        Fo3 = sf_output("out"); // wavetrans1
        Fo4 = sf_output("wavetrans2"); // wavetrans2 
        Fo5 = sf_output("amp1"); // amplitue1
        Fo6 = sf_output("amp2"); // amplitue2

		float **snap1, **snap2;
		float **sn1, **sn2;
		float *amp1, *amp2;

        float amax, amax1, amax2;

        /* Read/Write axes */
	    sf_axis az, ax;
	    az = sf_iaxa(Fo1,1); nz = sf_n(az); dz = sf_d(az)*1000;
	    ax = sf_iaxa(Fo1,2); nx = sf_n(ax); dx = sf_d(ax)*1000;
        fx=sf_o(ax);
        fz=sf_o(az);

        sf_warning("nx= %d nz= %d",nx,nz);
        sf_warning("dx= %f dz= %f",dx,dz);

		na = 720;
		da = 180.0/na;
		nb = nx*2;
		db = dx/2.0;

        sf_warning("da= %f db= %f",da,db);

		sf_putint(Fo3,"n1",nb);
		sf_putint(Fo3,"n2",na);
		sf_putfloat(Fo3,"d1",db);
		sf_putfloat(Fo3,"o1",0.0f);
		sf_putfloat(Fo3,"d2",da);
		sf_putfloat(Fo3,"o2",-90.0f);

		sf_putint(Fo4,"n1",nb);
		sf_putint(Fo4,"n2",na);
		sf_putfloat(Fo4,"d1",db);
		sf_putfloat(Fo4,"o1",0.0f);
		sf_putfloat(Fo4,"d2",da);
		sf_putfloat(Fo4,"o2",-90.0f);

		sf_putint(Fo5,"n2",1);
		sf_putint(Fo5,"n1",na);
		sf_putfloat(Fo5,"d1",da);
		sf_putfloat(Fo5,"o1",-90.0f);
		sf_putstring(Fo5,"label1","Phase angle");
		sf_putstring(Fo5,"unit1","Degree");

		sf_putint(Fo6,"n2",1);
		sf_putint(Fo6,"n1",na);
		sf_putfloat(Fo6,"d1",da);
		sf_putfloat(Fo6,"o1",-90.0f);
		sf_putstring(Fo6,"label1","Phase angle");
		sf_putstring(Fo6,"unit1","Degree");

		snap1=sf_floatalloc2(nz, nx);
		snap2=sf_floatalloc2(nz, nx);
		sn1=sf_floatalloc2(nb, na);
		sn2=sf_floatalloc2(nb, na);
		amp1=sf_floatalloc(na);
		amp2=sf_floatalloc(na);

		for(i=0;i<nx;i++){
            sf_floatread(snap1[i], nz, Fo1);
            sf_floatread(snap2[i], nz, Fo2);
		}
		for(i=0;i<na;i++){
		    amp1[i]=0.0;
		    amp2[i]=0.0;
			for(j=0;j<nb;j++){
				  sn1[i][j]=0.0;
				  sn2[i][j]=0.0;
			}
		}
		for(i=0;i<na;i++){
			float a=i*da;
			a *= SF_PI/180.0;
			for(j=0;j<nb;j++){
				float b=j*db;
				float x=(nx-1)*dx*0.5-b*cos(a);
				float z=(nz-1)*dz*0.5-b*sin(a);
				int ix, iz;
				ix=x/dx;
				iz=z/dz;
				if(ix>=0&&ix<nx&&iz>=0&&iz<nz/2+1){
			     //sf_warning("i=%d a=%f j=%d b=%f x=%f z=%f ix=%d iz=%d",i,a,j,b,x,z,ix,iz);
				  sn1[i][j]=snap1[ix][iz];
				  sn2[i][j]=snap2[ix][iz];
				}
			}
            sf_floatwrite(sn1[i], nb, Fo3);
            sf_floatwrite(sn2[i], nb, Fo4);
		}
		for(i=0;i<na;i++){
			amax=0.0;
		    for(j=0;j<nb;j++)
			   if(fabs(sn1[i][j])>amax) amax=fabs(sn1[i][j]);
			amp1[i]=amax;
            amax=0.0;
			for(j=0;j<nb;j++)
			   if(fabs(sn2[i][j])>amax) amax=fabs(sn2[i][j]);
			amp2[i]=amax;
        }
		
        amax1=0.0;
        amax2=0.0;
		for(i=na/2-2;i<na/2+2;i++){
			if(amp1[i]>amax1) amax1=amp1[i];
			if(amp2[i]>amax2) amax2=amp2[i];
		}
		for(i=0;i<na;i++){
            amp1[i] /= amax1;
            amp2[i] /= amax2;
		}

        sf_floatwrite(amp1, na, Fo5);
        sf_floatwrite(amp2, na, Fo6);
        free(*snap1);
        free(*snap2);
        free(*sn1);
        free(*sn2);
        free(amp1);
        free(amp2);
}
