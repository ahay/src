/* sort density*/
/*
  Copyright (C) 2010 Colorado School of Mines
  
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
#define EPS 3.

int main(int argc, char* argv[])
{
    bool verb;

    /* I/O files */
    sf_file Fin=NULL;
    sf_file Fout=NULL; 
    
  
    sf_axis az,ax,ay,aw;      /* cube axes */
    sf_axis dena1,dena2,dena3,dena4;
    float ****in;
    float **out;
   
    int nx,nz,ny,nw,nn;
    float dx,dz,dy,dw,ox,oz,oy,ow;
    float thrsh;
 

    int ii,jj,kk,ll;
    float x,y,z,w;
    int ind=0,n=0;
    float fx=0,fy=0,fz=0,fw=0;
    
    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getfloat("thrsh",&thrsh)) thrsh=0.05; /* verbosity flag */
    
    /*------------------------------------------------------------*/
    /* input files and axes */    
    Fin  = sf_input ("in");
   
    ax = sf_iaxa(Fin,1);
    ay = sf_iaxa(Fin,2);
    az = sf_iaxa(Fin,3);
    aw = sf_iaxa(Fin,4);
    nz=sf_n(az);
    nx=sf_n(ax);
    ny=sf_n(ay);
    nw=sf_n(aw);

    oz=sf_o(az);
    ox=sf_o(ax);
    oy=sf_o(ay);
    ow=sf_o(aw);

    dz=sf_d(az);
    dx=sf_d(ax);
    dy=sf_d(ay);
    dw=sf_d(aw);
    nn=nx*ny*nz*nw;

    in=sf_floatalloc4(nx,ny,nz,nw);   
    sf_floatread(in[0][0][0],nn,Fin); 
    if(verb) { 
	sf_raxa(ax); 
	sf_raxa(ay); 
	sf_raxa(az); 
	sf_raxa(aw); 
    }
    
   
    /*------------------------------------------------------------*/
   

    /*    /\*------------------------------------------------------------*\/ */
    /*     /\* output files*\/ */
    Fout = sf_output ("out");   
       
    for(ii=0;ii<nx;ii++){
	for(jj=0;jj<ny;jj++){
	    for(kk=0;kk<nz;kk++){
		for(ll=0;ll<nw;ll++){
		    if (ii<nx-1) fx=fabs(in[ll][kk][jj][ii+1]-in[ll][kk][jj][ii])/dx;
		    if (jj<ny-1) fy=fabs(in[ll][kk][jj+1][ii]-in[ll][kk][jj][ii])/dy;
		    if (kk<nz-1) fz=fabs(in[ll][kk+1][jj][ii]-in[ll][kk][jj][ii])/dz;
		    if (ll<nw-1) fw=fabs(in[ll+1][kk][jj][ii]-in[ll][kk][jj][ii])/dw;
		    
		    if(fx<EPS && fy<EPS && fz<EPS/90 && fw<EPS/180&& in[ll][kk][jj][ii]>thrsh  ){
			sf_warning("%f %f %f %f den=%f",
				   fx,fy,fz,fw,in[ll][kk][jj][ii] );		    
			x=ox+(ii+.5)*dx;
			y=oy+(jj+.5)*dy;
			z=oz+(kk+.5)*dz;
			w=ow+(ll+.5)*dw;
			sf_warning("%f %f %f %f den=%f",
				   x,y,z,w,in[ll][kk][jj][ii] );
			sf_warning("");
			n+=1;
		    }
		}
	    }
	}
    }
    out=sf_floatalloc2(n,4);
    dena1=sf_maxa( n, 0, 1);
    dena2=sf_maxa( 4, 0, 1);
    dena3=sf_maxa( 1, 0, 1);
    dena4=sf_maxa( 1, 0, 1);
    sf_oaxa(Fout,dena1,1);
    sf_oaxa(Fout,dena2,2);
    sf_oaxa(Fout,dena3,3);
    sf_oaxa(Fout,dena4,4);
    
    for(ii=0;ii<nx;ii++){
	for(jj=0;jj<ny;jj++){
	    for(kk=0;kk<nz;kk++){
		for(ll=0;ll<nw;ll++){
		    if (ii<nx-1) fx=fabs(in[ll][kk][jj][ii+1]-in[ll][kk][jj][ii])/dx;
		    if (jj<ny-1) fy=fabs(in[ll][kk][jj+1][ii]-in[ll][kk][jj][ii])/dy;
		    if (kk<nz-1) fz=fabs(in[ll][kk+1][jj][ii]-in[ll][kk][jj][ii])/dz;
		    if (ll<nw-1) fw=fabs(in[ll+1][kk][jj][ii]-in[ll][kk][jj][ii])/dw;
		    
		    if(fx<EPS && fy<EPS && fz<EPS && fw<EPS/180&& in[ll][kk][jj][ii]>thrsh  ){

/*		    if(in[ll][kk][jj][ii]>thrsh){*/
			out[0][ind]=ox+(ii+.5)*dx;
			out[1][ind]=oy+(jj+.5)*dy;
			out[2][ind]=oz+(kk+.5)*dz;
			out[3][ind]=ow+(ll+.5)*dw;
			
			ind+=1;
		
		    }
		}
	    }
	}
    }
     sf_floatwrite(out[0],n*4,Fout); 
     exit (0);
    
    
    
}


