/* 3D elastic time-domain FD modeling */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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
#ifdef _OPENMP
#include <omp.h>
#endif


int main(int argc, char* argv[])
{
    
    
    bool verb;
    int  nx,nz,nc,ic,iz,ix;
    int  ind=0;
    float tmp;
    float dx,dz,ox,oz;
    
    pt2d  *cc=NULL;
    float **eps=NULL,**del,***out=NULL;
    
    float x,y,lambda0,lambda1,lambda2,twoA;
    float lm00,lm01,lm02,lm10,lm11,lm12,lm20,lm21,lm22;
    
    float *a=NULL,*b=NULL,**llm=NULL;
    sf_axis ax,az,ac;
    sf_file Fin=NULL,Fdel=NULL, Fc=NULL,Fout=NULL; 
    
     /* init RSF */
    sf_init(argc,argv);
    

    if(! sf_getbool("verb",&verb)) verb=false;

    Fin = sf_input ("in" ); 
    Fdel= sf_input ("del" ); 
    Fout = sf_output ("out" );  
    Fc = sf_input ("cc" ); /* sample locations  */

    /* input axes */
    az = sf_iaxa(Fin,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
    ax = sf_iaxa(Fin,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);

    nz = sf_n(az); dz = sf_d(az); oz = sf_o(az);
    nx = sf_n(ax); dx = sf_d(ax); ox = sf_o(ax);

    eps = sf_floatalloc2(nz,nx);
    del = sf_floatalloc2(nz,nx);
    sf_floatread(eps[0],nz*nx,Fin);
    sf_floatread(del[0],nz*nx,Fdel);
     /* CIP coordinates */
    ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
    if(verb) sf_raxa(ac);
    nc = sf_n(ac);
    cc= (pt2d*) sf_alloc(nc,sizeof(*cc));
    pt2dread1(Fc,cc,nc,2); /* read coordinates */    
    /*  nc=3*/


     /* output files*/
    out = sf_floatalloc3(nz,nx,nc); 
    sf_oaxa(Fout,az,1);
    sf_oaxa(Fout,ax,2);
    sf_oaxa(Fout,ac,3);
    

    nc=3;
    a = sf_floatalloc(nc);
    b = sf_floatalloc(nc);
    llm = sf_floatalloc2(2,nc);

    /*find sample indexes*/
    for(ic=0; ic<nc; ic++) { 
	b[ic]=cc[ic].z; /*delta*/
	a[ic]=cc[ic].x; /*epsilon*/
	sf_warning("x%d=%f,z%d=%f: epsilon=%f delta=%f",ic,cc[ic].x,ic,cc[ic].z,a[ic],b[ic]);
    }
    
/*    a[0]=0.0; b[0]=0.0;*/
/*    a[2]=0.65; b[2]=0.0; */
/*    a[1]=0.0; b[1]=0.35;*/

    /*matrix*/
 
    twoA=(a[0]-a[2])*(b[1]-b[2])-(a[1]-a[2])*(b[0]-b[2]);
   /*  if(twoA<0) { */
/* 	tmp=a[1];a[1]=a[2];a[2]=tmp; */
/* 	tmp=b[1];b[1]=b[2];b[2]=tmp; */
/*     } */
   
    sf_warning("a=%f %f %f",a[0],a[1],a[2]);
    sf_warning("b=%f %f %f",b[0],b[1],b[2]);

    twoA=(a[0]-a[2])*(b[1]-b[2])-(a[1]-a[2])*(b[0]-b[2]);
    llm[0][0]=a[1]*b[2]-a[2]*b[1]; llm[0][1]=b[1]-b[2]; llm[0][2]=a[2]-a[1];
    llm[1][0]=a[2]*b[0]-a[0]*b[2]; llm[1][1]=b[2]-b[0]; llm[1][2]=a[0]-a[2];
    llm[2][0]=a[0]*b[1]-a[1]*b[0]; llm[2][1]=b[0]-b[1]; llm[2][2]=a[1]-a[0];
    sf_warning("");
    sf_warning("%f %f %f",llm[0][0], llm[0][1], llm[0][2]);
    sf_warning("%f %f %f",llm[1][0], llm[1][1], llm[1][2]);
    sf_warning("%f %f %f",llm[2][0], llm[2][1], llm[2][2]);
    sf_warning("area=%f",twoA);

    lm00=a[1]*b[2]-a[2]*b[1]; lm01=b[1]-b[2]; lm02=a[2]-a[1];
    lm10=a[2]*b[0]-a[0]*b[2]; lm11=b[2]-b[0]; lm12=a[0]-a[2];
    lm20=a[0]*b[1]-a[1]*b[0]; lm21=b[0]-b[1]; lm22=a[1]-a[0];
    sf_warning("");
    sf_warning("%f %f %f",lm00, lm01, lm02);
    sf_warning("%f %f %f",lm10, lm11, lm12);
    sf_warning("%f %f %f",lm20, lm21, lm22);








    /* find shape functions  */
    for    (ic=0; ic<nc; ic++) {
     for    (ix=0; ix<nx; ix++) {
	for(iz=0; iz<nz; iz++) {
	    out[ic][ix][iz]=0.;
	}
     }
    }

    if(twoA!=0){    
	for    (ix=0; ix<nx; ix++) {
	    for(iz=0; iz<nz; iz++) {
		x=eps[ix][iz];
		y=del[ix][iz];
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

}

