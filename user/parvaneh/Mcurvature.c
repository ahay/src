/* Curvature */

/*
  Copyright (C) 2004 University of Texas at Austin

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

int main (int argc, char *argv[])
{
    /*differentiation coefficient*/
    float c0=1./12., c1=1./6., c2=1./4., c3=1./9.;
    sf_file hor=NULL,cur=NULL;
    sf_axis ay,ax,at;
    int iy,ix,it;
    int ny,nx,nt;
    float idx,dx,dy,dt;
    const char *what;
    float ***a,***b,***c,***d,***e,***f,***val,***curv,***km,***kg,***kmax,***kmin,***kmp,***kmn,***kd,***ks,***kc;
    
    /*sf_init(!sf_getbool("verb",&verb)) verb=0;*/

    sf_init(argc,argv);
    hor=sf_input("in");
    cur=sf_output("out");

    at=sf_iaxa(hor,1); nt=sf_n(at); dt=sf_d(at);
    ay=sf_iaxa(hor,2); ny=sf_n(ay); dy=sf_d(ay);  
    ax=sf_iaxa(hor,3); nx=sf_n(ax); dx=sf_d(ax);

    sf_oaxa(cur,at,1);
    sf_oaxa(cur,ay,2);
    sf_oaxa(cur,ax,3);
   
    idx=1/dx;
    if (NULL == (what = sf_getstring("what"))) what="val";
    /* what to compute */
    val=sf_floatalloc3(nt,ny,nx); sf_floatread(val[0][0],nt*ny*nx,hor);

    a=sf_floatalloc3(nt,ny,nx);
    b=sf_floatalloc3(nt,ny,nx);
    c=sf_floatalloc3(nt,ny,nx);
    d=sf_floatalloc3(nt,ny,nx);
    e=sf_floatalloc3(nt,ny,nx);
    f=sf_floatalloc3(nt,ny,nx);
    curv=sf_floatalloc3(nt,ny,nx);
  
    switch (what[0]) {
	case 'm': /*Mean curvature*/
	    km=sf_floatalloc3(nt,ny,nx);
	    for (it=0; it<nt; it++){
		for (iy=1; iy<ny-1; iy++) {
		    for (ix=1; ix<nx-1; ix++) {
			a[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			b[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			c[ix][iy][it]=
			    idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			d[ix][iy][it]=
			    idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			e[ix][iy][it]=
			    idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			f[ix][iy][it]=
			    c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
			
			km[ix][iy][it]=
			    (a[ix][iy][it]*(1+e[ix][iy][it]*e[ix][iy][it])+b[ix][iy][it]*(1+d[ix][iy][it]*d[ix][iy][it])-(c[ix][iy][it]*d[ix][iy][it]*e[ix][iy][it]))/((1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*sqrtf(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));

			curv[ix][iy][it]=km[ix][iy][it];
		    }
		}
	    }
    
	    break;
        case 'g': /*Gaussian curvature*/  
	     kg=sf_floatalloc3(nt,ny,nx);
	     for (it=0; it<nt; it++){
		for (iy=1; iy<ny-1; iy++) {
		    for (ix=1; ix<nx-1; ix++) {
			a[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			b[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			c[ix][iy][it]=
			    idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			d[ix][iy][it]=
			    idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			e[ix][iy][it]=
			    idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			f[ix][iy][it]=
			    c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
	   
			kg[ix][iy][it]=
			    ((4*a[ix][iy][it]*b[ix][iy][it])-(c[ix][iy][it]*c[ix][iy][it]))/((1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			curv[ix][iy][it]=kg[ix][iy][it];
		    }
		}
	     }
	     break;
        case 'x': /*Maximum curvature*/
	    kmax=sf_floatalloc3(nt,ny,nx);
	    for (it=0; it<nt; it++){
		for (iy=1; iy<ny-1; iy++) {
		    for (ix=1; ix<nx-1; ix++) {
			a[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			b[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			c[ix][iy][it]=
			    idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			d[ix][iy][it]=
			    idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			e[ix][iy][it]=
			    idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			f[ix][iy][it]=
			    c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);

			km[ix][iy][it]=
			    (a[ix][iy][it]*(1+e[ix][iy][it]*e[ix][iy][it])+b[ix][iy][it]*(1+d[ix][iy][it]*d[ix][iy][it])-(c[ix][iy][it]*d[ix][iy][it]*e[ix][iy][it]))/((1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*sqrtf(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			
			kg[ix][iy][it]=
			    ((4*a[ix][iy][it]*b[ix][iy][it])-(c[ix][iy][it]*c[ix][iy][it]))/((1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			kmax[ix][iy][it]=
			    km[ix][iy][it]+sqrtf(km[ix][iy][it]*km[ix][iy][it]-kg[ix][iy][it]);
		    
			curv[ix][iy][it]=kmax[ix][iy][it];
		    }
		}
	    }
	    break;
	case 'i': /*Minimum curvature*/
	    kmin=sf_floatalloc3(nt,ny,nx);
	    for (it=0; it<nt; it++){
		for (iy=1; iy<ny-1; iy++) {
		    for (ix=1; ix<nx-1; ix++) {
			a[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			b[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			c[ix][iy][it]=
			    idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			d[ix][iy][it]=
			    idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			e[ix][iy][it]=
			    idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			f[ix][iy][it]=
			    c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
			km[ix][iy][it]=
			    (a[ix][iy][it]*(1+e[ix][iy][it]*e[ix][iy][it])+b[ix][iy][it]*(1+d[ix][iy][it]*d[ix][iy][it])-(c[ix][iy][it]*d[ix][iy][it]*e[ix][iy][it]))/((1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*sqrtf(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			
			kg[ix][iy][it]=
			    ((4*a[ix][iy][it]*b[ix][iy][it])-(c[ix][iy][it]*c[ix][iy][it]))/((1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));		   
			
			kmin[ix][iy][it]=
			    km[ix][iy][it]-sqrtf(km[ix][iy][it]*km[ix][iy][it]-kg[ix][iy][it]);
			curv[ix][iy][it]=kmin[ix][iy][it];
		    }
		}
	    }
	    break;
        case 'p': /*Most positive curvature*/
	    kmp=sf_floatalloc3(nt,ny,nx);
	    for (it=0; it<nt; it++){
		for (iy=1; iy<ny-1; iy++) {
		    for (ix=1; ix<nx-1; ix++) {
			a[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
			
			b[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			c[ix][iy][it]=
			    idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			d[ix][iy][it]=
			    idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			 
			e[ix][iy][it]=
			    idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);
			 
			f[ix][iy][it]=
			    c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
	    
	   
			kmp[ix][iy][it]=
			    (a[ix][iy][it]+b[ix][iy][it])-sqrtf((a[ix][iy][it]-b[ix][iy][it])*(a[ix][iy][it]-b[ix][iy][it])+(c[ix][iy][it]*c[ix][iy][it]));
			curv[ix][iy][it]=kmp[ix][iy][it];
		    }
		}
	    }
	    break;
        case 'n': /*Most negative curvature*/
	    kmn=sf_floatalloc3(nt,ny,nx);
	    for (it=0; it<nt; it++){
		for (iy=1; iy<ny-1; iy++) {
		    for (ix=1; ix<nx-1; ix++) {
			a[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			b[ix][iy][it]=
			    idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			c[ix][iy][it]=
			    idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			d[ix][iy][it]=
			    idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			e[ix][iy][it]=
			    idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);
			 
			f[ix][iy][it]=
			    c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
	   
			kmn[ix][iy][it]=
			    (a[ix][iy][it]+b[ix][iy][it])+sqrtf((a[ix][iy][it]-b[ix][iy][it])*(a[ix][iy][it]-b[ix][iy][it])+(c[ix][iy][it]*c[ix][iy][it]));
			curv[ix][iy][it]=kmn[ix][iy][it];
		    }
		}
	    }
	    break;
	 case 'd': /*Dip curvature*/
	     kd=sf_floatalloc3(nt,ny,nx);
	     for (it=0; it<nt; it++){
		 for (iy=1; iy<ny-1; iy++) {
		     for (ix=1; ix<nx-1; ix++) {
			 a[ix][iy][it]=
			     idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			 b[ix][iy][it]=
			     idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			 
			 c[ix][iy][it]=
			     idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			 d[ix][iy][it]=
			     idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			 e[ix][iy][it]=
			     idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			 f[ix][iy][it]=
			     c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
			
	    
			 kd[ix][iy][it]=
			     2*(a[ix][iy][it]*d[ix][iy][it]*d[ix][iy][it]+b[ix][iy][it]*e[ix][iy][it]*e[ix][iy][it]+c[ix][iy][it]*d[ix][iy][it]*e[ix][iy][it])/((d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*sqrtf(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			 curv[ix][iy][it]=kd[ix][iy][it];
		     }
		 }
	     }
	     break;
	  case 's': /*Strike curvature*/
	      ks=sf_floatalloc3(nt,ny,nx);
	      for (it=0; it<nt; it++){
		  for (iy=1; iy<ny-1; iy++) {
		      for (ix=1; ix<nx-1; ix++) {
			  a[ix][iy][it]=
			      idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			  b[ix][iy][it]=
			      idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			  c[ix][iy][it]=
			      idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
			  
			  d[ix][iy][it]=
			      idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			  e[ix][iy][it]=
			      idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			  f[ix][iy][it]=
			      c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
	   
	     
			  ks[ix][iy][it]=
			      2*(a[ix][iy][it]*e[ix][iy][it]*e[ix][iy][it]+b[ix][iy][it]*d[ix][iy][it]*d[ix][iy][it]-c[ix][iy][it]*d[ix][iy][it]*e[ix][iy][it])/((d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*sqrtf(1+d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			  curv[ix][iy][it]=ks[ix][iy][it];
		      }
		  }
	      }
	      break;
	 case 'c': /*Contour curvature*/
	     kc=sf_floatalloc3(nt,ny,nx);
	     for (it=0; it<nt; it++){
		 for (iy=1; iy<ny-1; iy++) {
		     for (ix=1; ix<nx-1; ix++) {
			 a[ix][iy][it]=
			     idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])-c1*(val[ix+1][iy][it]+val[ix][iy][it]+val[ix-1][iy][it]));
	    
			 b[ix][iy][it]=
			     idx*idx*(c0*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy][it]+val[ix-1][iy+1][it])-c1*(val[ix][iy-1][it]+val[ix][iy][it]+val[ix][iy+1][it]));
			   
			 c[ix][iy][it]=
			     idx*idx*c2*(val[ix+1][iy+1][it]+val[ix-1][iy-1][it]-val[ix+1][iy-1][it]-val[ix-1][iy+1][it]);
		   
			 d[ix][iy][it]=
			     idx*c1*(val[ix+1][iy+1][it]+val[ix][iy+1][it]+val[ix-1][iy+1][it]-val[ix+1][iy-1][it]-val[ix][iy-1][it]-val[ix-1][iy-1][it]);
			
			 e[ix][iy][it]=
			     idx*c1*(val[ix+1][iy-1][it]+val[ix+1][iy][it]+val[ix+1][iy+1][it]-val[ix-1][iy-1][it]-val[ix-1][iy][it]-val[ix-1][iy+1][it]);

			 f[ix][iy][it]=
			     c3*(2*(val[ix+1][iy][it]+val[ix][iy-1][it]+val[ix][iy+1][it]+val[ix-1][iy][it])-(val[ix+1][iy-1][it]+val[ix+1][iy+1][it]+val[ix-1][iy-1][it]+val[ix-1][iy+1][it])+5*val[ix][iy][it]);
	    
			 kc[ix][iy][it]=
			     2*(a[ix][iy][it]*e[ix][iy][it]*e[ix][iy][it]+b[ix][iy][it]*d[ix][iy][it]*d[ix][iy][it]-c[ix][iy][it]*d[ix][iy][it]*e[ix][iy][it])/((d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it])*sqrt(d[ix][iy][it]*d[ix][iy][it]+e[ix][iy][it]*e[ix][iy][it]));
			 curv[ix][iy][it]=kc[ix][iy][it];
		     }
		 }
	     }
	     break;
   
		     
    }

    sf_floatwrite(curv[0][0],nt*ny*nx,cur);

    sf_close();

    exit (0);
}    

    
		 
	       
