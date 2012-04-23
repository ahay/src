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

/* Roberts, A., 2001, Curvature attributes and their application to
   3-D interpreted horizons: First Break, 19, no. 2, 85–99. */

int main (int argc, char *argv[])
{
    /*differentiation coefficient*/
    float c0=1./12., c1=1./6., c2=1./4.;
    sf_file hor=NULL,cur=NULL;
    sf_axis ay,ax,at;
    int iy,ix,it;
    int ny,nx,nt;
    float idx,dx,idy,dy;
    const char *what;
    float a,b,c,d,e,***val,km,kg,**slice;
    
    sf_init(argc,argv);
    hor=sf_input("in");
    cur=sf_output("out");

    at=sf_iaxa(hor,1); nt=sf_n(at);
    ay=sf_iaxa(hor,2); ny=sf_n(ay); dy=sf_d(ay);  
    ax=sf_iaxa(hor,3); nx=sf_n(ax); dx=sf_d(ax);
   
    idx=1.0/dx;
    idy=1.0/dy;

    if (NULL == (what = sf_getstring("what"))) what="val";
    /* what to compute */
    val =sf_floatalloc3(nt,ny,nx);
    slice = sf_floatalloc2(ny,nx);

    sf_floatread(val[0][0],nt*ny*nx,hor);
  
    for (it=0; it<nt; it++){
	for (ix=0; ix<nx; ix++) {
	    for (iy=0; iy<ny; iy++) {
		slice[ix][iy] = val[ix][iy][it];
		val[ix][iy][it] = 0.0;
	    }
	}


	for (ix=1; ix<nx-1; ix++) {
	    for (iy=1; iy<ny-1; iy++) {
		a=idy*idy*(c0*(slice[ix+1][iy-1]+slice[ix+1][iy+1]+slice[ix][iy-1]+slice[ix][iy+1]+slice[ix-1][iy-1]+slice[ix-1][iy+1])-c1*(slice[ix+1][iy]+slice[ix][iy]+slice[ix-1][iy]));
	    
		b=idx*idx*(c0*(slice[ix+1][iy-1]+slice[ix+1][iy]+slice[ix+1][iy+1]+slice[ix-1][iy-1]+slice[ix-1][iy]+slice[ix-1][iy+1])-c1*(slice[ix][iy-1]+slice[ix][iy]+slice[ix][iy+1]));
		
		c=idx*idy*c2*(slice[ix+1][iy+1]+slice[ix-1][iy-1]-slice[ix+1][iy-1]-slice[ix-1][iy+1]);
		   
		d=idy*c1*(slice[ix+1][iy+1]+slice[ix][iy+1]+slice[ix-1][iy+1]-slice[ix+1][iy-1]-slice[ix][iy-1]-slice[ix-1][iy-1]);
			
		e=idx*c1*(slice[ix+1][iy-1]+slice[ix+1][iy]+slice[ix+1][iy+1]-slice[ix-1][iy-1]-slice[ix-1][iy]-slice[ix-1][iy+1]);
			
		switch (what[0]) {
		    case 'm': /*Mean curvature*/
			val[ix][iy][it]=(a*(1+e*e)+b*(1+d*d)-(c*d*e))/((1+d*d+e*e)*sqrtf(1+d*d+e*e));
			break;
		    case 'g': /*Gaussian curvature*/  
			val[ix][iy][it]=((4*a*b)-(c*c))/((1+d*d+e*e)*(1+d*d+e*e));
			break;
		    case 'x': /*Maximum curvature*/
			km=(a*(1+e*e)+b*(1+d*d)-(c*d*e))/((1+d*d+e*e)*sqrtf(1+d*d+e*e));		      
			kg=((4*a*b)-(c*c))/((1+d*d+e*e)*(1+d*d+e*e));
			val[ix][iy][it]=km+sqrtf(km*km-kg);
			break;
		    case 'i': /*Minimum curvature*/
			km=(a*(1+e*e)+b*(1+d*d)-(c*d*e))/((1+d*d+e*e)*sqrtf(1+d*d+e*e));		      
			kg=((4*a*b)-(c*c))/((1+d*d+e*e)*(1+d*d+e*e));		   			
			val[ix][iy][it]=km-sqrtf(km*km-kg);
			break;
		    case 'p': /*Most positive curvature*/
			val[ix][iy][it]=(a+b)-sqrtf((a-b)*(a-b)+(c*c));
			break;
		    case 'n': /*Most negative curvature*/
			val[ix][iy][it]=(a+b)+sqrtf((a-b)*(a-b)+(c*c));		
			break;
		    case 'd': /*Dip curvature*/
			val[ix][iy][it]=2*(a*d*d+b*e*e+c*d*e)/((d*d+e*e)*(1+d*d+e*e)*sqrtf(1+d*d+e*e));
			break;
		    case 's': /*Strike curvature*/
			val[ix][iy][it]=2*(a*e*e+b*d*d-c*d*e)/((d*d+e*e)*sqrtf(1+d*d+e*e));		      
			break;
		    case 'c': /*Contour curvature*/
			val[ix][iy][it]=2*(a*e*e+b*d*d-c*d*e)/((d*d+e*e)*sqrt(d*d+e*e));
			break;
		}
	    }
	}
    }

    sf_floatwrite(val[0][0],nt*ny*nx,cur);
    exit (0);
}

