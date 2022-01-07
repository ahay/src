/* */
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

#include <math.h>
#define INF(x,y,z) x+y+z

static int fx(int x,int ix,int iy,int iz,int ***g)
{
    return (x-ix)*(x-ix)+g[iy][ix][iz]*g[iy][ix][iz];
}


static int fy(int y,int iy,int ix,int iz,int ***g)
{
    return (y-iy)*(y-iy)+g[iy][ix][iz];
}


static int Sepx(int i, int u, int ***g, int iy, int iz)
{
    float tmp;
    tmp=(u*u-i*i+
	 g[iy][u][iz]*g[iy][u][iz]-
	 g[iy][i][iz]*g[iy][i][iz]) / (2*(u-i));    
    return (int)tmp;
}


static int Sepy(int i, int u, int ***g, int ix, int iz)
{
    float tmp;
    tmp=(u*u-i*i+
	 g[u][ix][iz]-
	 g[i][ix][iz]) / (2*(u-i));    
    return (int)tmp;
}


int main(int argc, char* argv[])
{
    bool verb;

    /* I/O files */
    sf_file Fin=NULL; 
    sf_file Fout=NULL; 
    
  
    sf_axis az,ax,ay;      /* cube axes */
    float ***in;              
    float ***out;
    int ***g;
    int ***h;
    int ix,iy,iz,nx,nz,ny,nn;

    int q,u,w;
    int *s,*t;
     
    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    
    
    /*------------------------------------------------------------*/
    /* input files and axes */    
    Fin = sf_input ("in");

    az = sf_iaxa(Fin,1);
    ax = sf_iaxa(Fin,2);
    ay = sf_iaxa(Fin,3);
    nz=sf_n(az);
    nx=sf_n(ax);
    ny=sf_n(ay);

    nn=sf_n(ax)*sf_n(az)*sf_n(ay);

    in=sf_floatalloc3(nz,nx,ny);   
    sf_floatread(in[0][0],nn,Fin); 
    if(verb) { 
	sf_raxa(az); 
	sf_raxa(ax); 
	sf_raxa(ay); 
    }
    
    /*    /\*------------------------------------------------------------*\/ */
    /*     /\* output files*\/ */
    Fout = sf_output ("out");
    
    out=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    g=sf_intalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    h=sf_intalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    
    sf_oaxa(Fout,az,1);
    sf_oaxa(Fout,ax,2);
    sf_oaxa(Fout,ay,3);
    /*------------------------------------------------------------*/
    /* scan axis 1 */
    for (iy=0;iy<ny;iy++){		
	for (ix=0;ix<nx;ix++){
	    /* scan 1    */		
	    if (in[iy][ix][0]) 
		g[iy][ix][0]=0;
	    else 
		g[iy][ix][0]=INF(nx,ny,nz);
	    for (iz=1;iz<nz;iz++){
		if (in[iy][ix][iz]) 
		    g[iy][ix][iz]=0;
		else 
		    g[iy][ix][iz]=g[iy][ix][iz-1]+1;
	    }
	      /* scan 2    */
	    for (iz=nz-2;iz>=0;iz--){
		if (g[iy][ix][iz+1]<g[iy][ix][iz])
		    g[iy][ix][iz]=g[iy][ix][iz+1]+1;
	    }

	}
    }



    /* scan axis 2 */
    s=sf_intalloc(nx);
    t=sf_intalloc(nx);   
   
    for (iy=0;iy<ny;iy++){
	for (iz=0;iz<nz;iz++){	
	    q=0;s[0]=0;t[0]=0;
	    /* scan 3*/
	    for (u=1;u<nx;u++){
		while (q>=0 && fx(t[q],s[q],iy,iz,g) > fx(t[q],u,iy,iz,g) )
		{q=q-1;}
		if (q<0) {
		    q=0;s[0]=u;
		}
		else{
		    w=Sepx(s[q],u,g,iy,iz)+1;
		    if (w<nx){
			q=q+1;
			s[q]=u;
			t[q]=w;
		    }
		}
	    }	    
	    /* scan 4*/
	    for (u=nx-1;u>=0;u--){
		h[iy][u][iz]=fx(u,s[q],iy,iz,g);
		if (u==t[q]) q=q-1;		    
	    }
	    
	}
    }


    /* scan axis 3 */
    s=sf_intalloc(ny);
    t=sf_intalloc(ny);

    for (ix=0;ix<nx;ix++){
	for (iz=0;iz<nz;iz++){
	    q=0;s[0]=0;t[0]=0;
	    /* scan 3*/
	    for (u=1;u<ny;u++){
		while (q>=0 && fy(t[q],s[q],ix,iz,h) > fy(t[q],u,ix,iz,h) )
		{q=q-1;}
		if (q<0) {
		    q=0;s[0]=u;
		}
		else{
		    w=Sepy(s[q],u,h,ix,iz)+1;
		    if (w<ny){
			q=q+1;
			s[q]=u;
			t[q]=w;
		    }
		}
	    }
	    /* scan 4*/
	    for (u=ny-1;u>=0;u--){
		out[u][ix][iz]=fy(u,s[q],ix,iz,h);
		if (u==t[q]) q=q-1;
	    }
	    
	}
    }




  
    sf_floatwrite(out[0][0],nn,Fout); 
    exit (0);
    
    
    
}


