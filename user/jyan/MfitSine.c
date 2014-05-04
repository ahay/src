/* fit a shifted sine curve*/
/*
  Copyright (C) 2009 Colorado School of Mines
  
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

#include <math.h>

int main(int argc, char* argv[])
{
    
    
    bool verb;
    int  nx,nz,nc,iz,ix;
    float nu1,nu2,nu3;

    
    float **nu=NULL,***out=NULL;
    float *nuref=NULL;
    sf_axis ax,az,ac;
    sf_file Fin=NULL,Fout=NULL; 
    
     /* init RSF */
    sf_init(argc,argv);
    

    if(! sf_getbool("verb",&verb)) verb=false;

    if(! sf_getfloat("nu1",&nu1)) nu1=0;
    if(! sf_getfloat("nu2",&nu2)) nu2=SF_PI/2; 
    if(! sf_getfloat("nu3",&nu3)) nu3=-SF_PI/2; 

    sf_warning("nu1=%f nu2=%f n3=%f",nu1,nu2,nu3);

    Fin = sf_input ("in" ); 
    Fout = sf_output ("out" );  

    nc=3;
    /* input axes */
    az = sf_iaxa(Fin,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
    ax = sf_iaxa(Fin,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);
    ac = sf_maxa(nc,0,1);
	
    nz = sf_n(az); /* dz = sf_d(az); oz = sf_o(az); */
    nx = sf_n(ax); /* dx = sf_d(ax); ox = sf_o(ax); */

    nu = sf_floatalloc2(nz,nx); 
    sf_floatread(nu[0],nz*nx,Fin);

    /* CIP coordinates */

    /*  nc=3*/

     /* output files*/
    
    out = sf_floatalloc3(nz,nx,nc);
    sf_oaxa(Fout,az,1);
    sf_oaxa(Fout,ax,2);
    sf_oaxa(Fout,ac,3);
    
    float *a=NULL,*b=NULL;

    a = sf_floatalloc(nc);
    b = sf_floatalloc(nc);
    nuref = sf_floatalloc(nc);
    nuref[0]=nu1;
    nuref[1]=nu2;
    nuref[2]=nu3;
    int ic;
    for(ic=0; ic<nc; ic++) { 
	b[ic]=sin(2*nuref[ic]*SF_PI/180); /*delta*/
	a[ic]=cos(2*nuref[ic]*SF_PI/180); /*epsilon*/
	sf_warning("x=%f y=%f",a[ic],b[ic]);
    }
    float x,y,lambda0,lambda1,lambda2,twoA;
    float lm00,lm01,lm02,lm10,lm11,lm12,lm20,lm21,lm22;
    twoA=(a[0]-a[2])*(b[1]-b[2])-(a[1]-a[2])*(b[0]-b[2]);
    lm00=a[1]*b[2]-a[2]*b[1]; lm01=b[1]-b[2]; lm02=a[2]-a[1];
    lm10=a[2]*b[0]-a[0]*b[2]; lm11=b[2]-b[0]; lm12=a[0]-a[2];
    lm20=a[0]*b[1]-a[1]*b[0]; lm21=b[0]-b[1]; lm22=a[1]-a[0];
    sf_warning("");
    sf_warning("%f %f %f",lm00, lm01, lm02);
    sf_warning("%f %f %f",lm10, lm11, lm12);
    sf_warning("%f %f %f",lm20, lm21, lm22);



    for    (ic=0; ic<nc; ic++) {
	for    (ix=0; ix<nx; ix++) {
	    for(iz=0; iz<nz; iz++) {
		out[ic][ix][iz]=0.;
	    }
	}
    }
    sf_warning("nc=%d",nc);
    if(twoA!=0){    
	for    (ix=0; ix<nx; ix++) {
	    for(iz=0; iz<nz; iz++) {
		x=cos(2*nu[ix][iz]*SF_PI/180);
		y=sin(2*nu[ix][iz]*SF_PI/180);
/*		sf_warning("x,y=%f %f",x,y);*/
		lambda0=lm00 + lm01*x + lm02*y;
		lambda1=lm10 + lm11*x + lm12*y;
		lambda2=lm20 + lm21*x + lm22*y;
		out[0][ix][iz]=lambda0/twoA;
		out[1][ix][iz]=lambda1/twoA;
		out[2][ix][iz]=lambda2/twoA;
	    }
	}
    }
    sf_floatwrite(out[0][0],nz*nx*nc,Fout);


    exit (0);
    
    
    
}
