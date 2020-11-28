/* Curvature in stratigraphic coordinates */

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

#include "curvature.h"

/* Roberts, A., 2001, Curvature attributes and their application to
   3-D interpreted horizons: First Break, 19, no. 2, 85–99. */

int main (int argc, char *argv[])
{
    /*differentiation coefficient*/    
    sf_file xhor=NULL,cur=NULL, yhor=NULL, zhor=NULL;
    sf_axis ay,ax,at;
    int iy,ix,it;
    int ny,nx,nt;
    float idx,dx,idy,dy, scale;
    const char *what;    
    float xuu,yuu,zuu,xvv,yvv,zvv,***valx, ***valy, ***valz,km,kg,**slicex, **slicey, **slicez;
    float xu,yu,zu,xv,yv,zv,n1,n2,n3,nn1,nn2,nn3,ee,ff,gg,e,f,g;
    float zuv, xuv, yuv;

    sf_init(argc,argv);
    xhor=sf_input("in");
    yhor=sf_input("yh");
    zhor=sf_input("zh");
    cur=sf_output("out");

    at=sf_iaxa(xhor,1); nt=sf_n(at); 
    ay=sf_iaxa(xhor,2); ny=sf_n(ay); dy=sf_d(ay);  
    ax=sf_iaxa(xhor,3); nx=sf_n(ax); dx=sf_d(ax);
   
    if (NULL == (what = sf_getstring("what"))) what="valx";
    /* what to compute */

    
    if (!sf_getfloat("scale",&scale)) scale=1.0;
    /* scaling (from time to depth) */

   

    idx=1.0/dx;
    idy=1.0/dy;

    valx =sf_floatalloc3(nt,ny,nx);
    valy =sf_floatalloc3(nt,ny,nx);
    valz =sf_floatalloc3(nt,ny,nx);
    slicex = sf_floatalloc2(ny,nx);
    slicey = sf_floatalloc2(ny,nx);
    slicez = sf_floatalloc2(ny,nx);
    
    
    sf_floatread(valx[0][0],nt*ny*nx,xhor);
    sf_floatread(valy[0][0],nt*ny*nx,yhor);
    sf_floatread(valz[0][0],nt*ny*nx,zhor);

    for (it=0; it<nt; it++){
	for (ix=0; ix<nx; ix++) {
	    for (iy=0; iy<ny; iy++) {
		slicex[ix][iy] = valx[ix][iy][it];
		valx[ix][iy][it] = 0.0;

		slicey[ix][iy] = valy[ix][iy][it];
		valy[ix][iy][it] = 0.0;

		slicez[ix][iy] = scale*valz[ix][iy][it];
		valz[ix][iy][it] = 0.0;
	    }
	}
	
	
	for (ix=1; ix<nx-1; ix++) {
	    for (iy=1; iy<ny-1; iy++) {
		xuu=idy*idy*fd1_2(slicex, iy, ix, 'U');
		yuu=idy*idy*fd1_2(slicey, iy, ix, 'U');
		zuu=idy*idy*fd1_2(slicez, iy, ix, 'U');

		xvv=idx*idx*fd1_2(slicex, iy, ix, 'V');
		yvv=idx*idx*fd1_2(slicey, iy, ix, 'V');
		zvv=idx*idx*fd1_2(slicez, iy, ix, 'V');
		
		xuv=idy*idx*fd1_2(slicex, iy, ix, 'W');
		yuv=idy*idx*fd1_2(slicey, iy, ix, 'W');
		zuv=idy*idx*fd1_2(slicez, iy, ix, 'W');		
	
		xu=idy*fd1_2(slicex, iy, ix, 'u');
		yu=idy*fd1_2(slicey, iy, ix, 'u');
		zu=idy*fd1_2(slicez, iy, ix, 'u');

		xv=idx*fd1_2(slicex, iy, ix, 'v');
		yv=idx*fd1_2(slicey, iy, ix, 'v');
		zv=idx*fd1_2(slicez, iy, ix, 'v');
		
		ee=(xu*xu)+(yu*yu)+(zu*zu);
		ff=(xu*xv)+(yu*yv)+(zu*zv);
		gg=(xv*xv)+(yv*yv)+(zv*zv);
		
	       	n1=(yu*zv)-(zu*yv);
		n2=(zu*xv)-(xu*zv);
		n3=(xu*yv)-(yu*xv);
		
		nn1=n1/(sqrtf(n1*n1+n2*n2+n3*n3)+1E-20);
		nn2=n2/(sqrtf(n1*n1+n2*n2+n3*n3)+1E-20);
		nn3=n3/(sqrtf(n1*n1+n2*n2+n3*n3)+1E-20);
		
		e=(xuu*nn1)+(yuu*nn2)+(zuu*nn3);
		f=(xuv*nn1)+(yuv*nn2)+(zuv*nn3);
		g=(xvv*nn1)+(yvv*nn2)+(zvv*nn3);

		switch (what[0]) {
		    case 'm': /*Mean curvature*/
			valx[ix][iy][it]=(e*gg-2*f*ff+g*ee)/(2*(ee*gg-ff*ff)+1E-20);
			break;
		    case 'g': /*Gaussian curvature*/  
			valx[ix][iy][it]=(e*g-f*f)/(ee*gg-ff*ff+1E-20); /* flipped sign */
			break;
		    case 'x': /*Maximum curvature*/
			km=(e*gg-2*f*ff+g*ee)/(2*(ee*gg-ff*ff)+1E-20);		      
			kg=(e*g-f*f)/(ee*gg-ff*ff+1E-20);
			valx[ix][iy][it]=km+sqrtf(km*km-kg);
			break;
		    case 'i': /*Minimum curvature*/
			km=(e*gg-2*f*ff+g*ee)/(2*(ee*gg-ff*ff)+1E-20);		      
			kg=(e*g-f*f)/(ee*gg-ff*ff+1E-20);		   			
			valx[ix][iy][it]=-km+sqrtf(km*km-kg); /* flipped sign */
			break;		   
			
		}
	    }
	}
    }
    
        
    
    sf_floatwrite(valx[0][0],nt*ny*nx,cur);
    exit (0);
}

