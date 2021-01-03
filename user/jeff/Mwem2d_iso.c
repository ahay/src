/* 2D ISOTROPIC wave-equation finite-difference migration */
/*
  Copyright (C) 2012 The University of Western Australia
  
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

#include "wem2d_iso_wemig.h"

int main(int argc, char* argv[])
{
    int ns,nw,nx,nz,nh,is;
    float ox,oz,dx,dz,ow,dw,oh,dh;
    bool adj,add,wantwf;
    /* cube axes */
    sf_axis ax,aw,az,as,ah,anull;
    int nxtap;
    bool verbose;
    /* I/O arrays */
    float **vel=NULL;  /* v array */
    float ***xig=NULL; /* output XIG array */
    sf_complex **rwf=NULL; /* receiver wavefield array */
    sf_complex **swf=NULL; /* Source wavefield array */
  
    /* I/O files */
    sf_file Fvel  = NULL; /* velocity file */
    sf_file Fxig  = NULL; /* input XIG file */
    sf_file Fswf  = NULL; /* input SWF at iz=0 */
    sf_file Frwf  = NULL; /* input RWF at iz=0 */
    sf_file Fxigo = NULL; /* output xig file */
    sf_file Fswfo = NULL; /* output SWF at iz=NZ-1 file*/
    sf_file Frwfo = NULL; /* output RWF at iz=NZ-1 file*/

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    
    /* Set up file objects */
    Fvel = sf_input("vel"); /* input velocity */

    Fxig = sf_input("in"); /* input xig */
    Fswf = sf_input("swf"); sf_settype(Fswf,SF_COMPLEX);/* INPUT SWF at iz=0 */
    Frwf = sf_input("rwf"); sf_settype(Frwf,SF_COMPLEX);/* INPUT RWF at iz=0 */
 
    Fxigo = sf_output("out"); /* OUTPUT XIG */
    Fswfo = sf_output("swfout"); sf_settype(Fswfo,SF_COMPLEX);/* OUTPUT SWF at iz=NZ-1 */
    Frwfo = sf_output("rwfout"); sf_settype(Frwfo,SF_COMPLEX);/* OUTPUT rWF at iz=NZ-1 */

    /* Read in parameter */ 
    if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
    if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */
    if(! sf_getbool("adj",&adj)) adj=true; /* ADJOINT flag */
    sf_warning("Performing ADJ=%s",(adj)?"true":"false");
    if(! sf_getbool("add",&add)) add=false;
    if(! sf_getbool("wantwf",&wantwf)) wantwf=false; /* Want output wavefields */

    /* Get axes information */
    ax = sf_iaxa(Frwf,1);
    aw = sf_iaxa(Frwf,2);
    as = sf_iaxa(Frwf,3);
    az = sf_iaxa(Fvel,2);
    ah = sf_iaxa(Frwf,1);
    anull=sf_iaxa(Frwf,1);

    /* Pull out individual values */
    nx = sf_n(ax); ox = sf_o(ax); dx = sf_d(ax);
    nw = sf_n(aw); ow = sf_o(aw); dw = sf_d(aw);
    ns = sf_n(as);
    nz = sf_n(az); oz = sf_o(az); dz = sf_d(az);

    /* xig shift parameter */
    if(! sf_getint("nh",&nh)) nh=0;
    nh = 2*nh+1;
    dh = dx;
    oh = -(float)nh*dh+dx;

    /* Set up xig shift axis */
    sf_setn(ah,nh);
    sf_seto(ah,oh);
    sf_setd(ah,dh);
    sf_setn(anull,1);
    sf_setd(anull,1.);
    sf_seto(anull,0.);

    /* Output axes information */
    if (wantwf) {
	sf_oaxa(Fswfo,ax,1);
	sf_oaxa(Fswfo,aw,2);
	sf_oaxa(Fswfo,anull,3);

	sf_oaxa(Frwfo,ax,1);
	sf_oaxa(Frwfo,aw,2);
	sf_oaxa(Frwfo,anull,3);
    }
 
    /* Allocation */
    vel =   sf_floatalloc2(nx,nz);
    xig =   sf_floatalloc3(nh,nx,nz);
    rwf = sf_complexalloc2(nx,nw);
    swf = sf_complexalloc2(nx,nw);
  
    /* Read in velocity */
    sf_floatread( vel[0],nx*nz,Fvel );	
  
    /* Read in xig */
    sf_floatread( xig[0][0],nx*nz*nh,Fxig );

    /* Initialize propagation module */
    sf_raxa(az);
    sf_raxa(ax);
    sf_raxa(aw);
    sf_raxa(ah);
  
    wem2d_iso_wemig_init(nx,nz,nw,nh,ox,oz,ow,oh,dx,dz,dw,dh,nxtap,vel);
  
    /* Loop over inline shots */
    for (is=0; is < ns; is++) {
    
	/* Read in SWF and RWF  */
	sf_complexread(swf[0],nx*nw,Fswf);
	sf_complexread(rwf[0],nx*nw,Frwf);	
    
	/* IMAGE INDIVIDUAL SHOTS */
	wem2d_iso_wemig(adj,add,vel,swf,rwf,xig);	

	if (wantwf) {
	    /* Write out SWF and RWF */
	    sf_warning("WRITING OUT RWF and SWF");
	    sf_complexwrite(swf[0],nx*nw,Fswfo);
	    sf_complexwrite(rwf[0],nx*nw,Frwfo);
	}
    }
  
    sf_warning("WRITING OUT XIG");
    /* WRITE OUT XIG */
    sf_floatwrite(xig[0][0],nh*nx*nz,Fxigo);

    wem2d_iso_close();
  
    exit(0);
}
