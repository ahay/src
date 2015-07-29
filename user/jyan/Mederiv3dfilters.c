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
#include "tool3PSS.h"
#include "dsyevv3.h"
#include "dsyevc3.h"


#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;
    bool stat;

    char *domain;
    
    /* I/O files */
    sf_file Fccc=NULL;
    sf_file Fspk=NULL;  /* template for derivative operator */    
    sf_file Fpzdel=NULL; /* pz derivative */
    sf_file Fpxdel=NULL; /* px derivative */
    sf_file Fpydel=NULL; /* py derivative */

    sf_file Fvzdel=NULL; /* vz derivative */
    sf_file Fvxdel=NULL; /* vx derivative */
    sf_file Fvydel=NULL; /* vy derivative */
    
    
    sf_file Fhzdel=NULL; /* hz derivative */
    sf_file Fhxdel=NULL; /* hx derivative */
    sf_file Fhydel=NULL; /* hy derivative */

    sf_axis ax,ay,az;      /* cube axes */
    sf_axis sx,sy,sz;      /* cube axes */
    
    float ***c11,***c12,***c13,***c22,***c23,***c33,***c44,***c55,***c66; /* input file */
    /* float ***tmp; */
/*    float  mc11, mc12, mc13, mc22, mc23, mc33, mc44, mc55, mc66;*/
    float ***pxdel,***pydel,***pzdel;                /* output derivative */
    float ***vxdel,***vydel,***vzdel;                /* output derivative */
    float ***hxdel,***hydel,***hzdel;                /* output derivative */
    int ix,iy,iz,nnn,ndel;
    wfs3d wfs3;
    
    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("stat",&stat)) stat=true;  /* stationary operator */
    

    
    if(NULL== (domain=sf_getstring("domain"))) domain="k";
    if(verb) sf_warning("domain= %s",domain);
    /*------------------------------------------------------------*/
    /* input files and axes */
    Fccc = sf_input ("ccc");      
    az = sf_iaxa(Fccc,1);
    ax = sf_iaxa(Fccc,2);
    ay = sf_iaxa(Fccc,3);


    Fspk = sf_input ("in");
    sz = sf_iaxa(Fspk,1);
    sx = sf_iaxa(Fspk,2);
    sy = sf_iaxa(Fspk,3);

    if(verb) { 
	sf_raxa(sz); 
	sf_raxa(sx); 
	sf_raxa(sy); 
	sf_raxa(az); 
	sf_raxa(ax); 
	sf_raxa(ay);
    }

    /*------------------------------------------------------------*/
    /* output files*/
    Fpzdel = sf_output ("pzdel");
    Fpxdel = sf_output ("pxdel");
    Fpydel = sf_output ("pydel");

    Fhzdel = sf_output ("hzdel");
    Fhxdel = sf_output ("hxdel");
    Fhydel = sf_output ("hydel");
    
    Fvzdel = sf_output ("vzdel");
    Fvxdel = sf_output ("vxdel");
    Fvydel = sf_output ("vydel");

    if(domain[0]=='k'){
	sf_oaxa(Fpxdel,az,1);sf_oaxa(Fpxdel,ax,2);sf_oaxa(Fpxdel,ay,3);
	sf_oaxa(Fpydel,az,1);sf_oaxa(Fpydel,ax,2);sf_oaxa(Fpydel,ay,3);	
	sf_oaxa(Fpzdel,az,1);sf_oaxa(Fpzdel,ax,2);sf_oaxa(Fpzdel,ay,3);
	pzdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	pydel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	pxdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));

	sf_oaxa(Fvxdel,az,1);sf_oaxa(Fvxdel,ax,2);sf_oaxa(Fvxdel,ay,3);
	sf_oaxa(Fvydel,az,1);sf_oaxa(Fvydel,ax,2);sf_oaxa(Fvydel,ay,3);	
	sf_oaxa(Fvzdel,az,1);sf_oaxa(Fvzdel,ax,2);sf_oaxa(Fvzdel,ay,3);
	vzdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	vydel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	vxdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));

	sf_oaxa(Fhxdel,az,1);sf_oaxa(Fhxdel,ax,2);sf_oaxa(Fhxdel,ay,3);
	sf_oaxa(Fhydel,az,1);sf_oaxa(Fhydel,ax,2);sf_oaxa(Fhydel,ay,3);	
	sf_oaxa(Fhzdel,az,1);sf_oaxa(Fhzdel,ax,2);sf_oaxa(Fhzdel,ay,3);
	hzdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	hydel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	hxdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));

	wfs3 = wfsep3_init(ax,ay,az);
    }
    else{
	sf_oaxa(Fpxdel,sz,1);sf_oaxa(Fpxdel,sx,2);sf_oaxa(Fpxdel,sy,3);
	sf_oaxa(Fpydel,sz,1);sf_oaxa(Fpydel,sx,2);sf_oaxa(Fpydel,sy,3);
	sf_oaxa(Fpzdel,sz,1);sf_oaxa(Fpzdel,sx,2);sf_oaxa(Fpzdel,sy,3);
	pzdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	pydel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	pxdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));

	sf_oaxa(Fhxdel,sz,1);sf_oaxa(Fhxdel,sx,2);sf_oaxa(Fhxdel,sy,3);
	sf_oaxa(Fhydel,sz,1);sf_oaxa(Fhydel,sx,2);sf_oaxa(Fhydel,sy,3);
	sf_oaxa(Fhzdel,sz,1);sf_oaxa(Fhzdel,sx,2);sf_oaxa(Fhzdel,sy,3);
	hzdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	hydel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	hxdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	
	sf_oaxa(Fvxdel,sz,1);sf_oaxa(Fvxdel,sx,2);sf_oaxa(Fvxdel,sy,3);
	sf_oaxa(Fvydel,sz,1);sf_oaxa(Fvydel,sx,2);sf_oaxa(Fvydel,sy,3);
	sf_oaxa(Fvzdel,sz,1);sf_oaxa(Fvzdel,sx,2);sf_oaxa(Fvzdel,sy,3);
	vzdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	vydel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	vxdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	

	wfs3 = wfsep3_init(sx,sy,sz);
    }
        
    if(!stat) {
	sf_oaxa(Fpzdel,az,4);sf_oaxa(Fpzdel,ax,5);sf_oaxa(Fpzdel,ay,6);
	sf_oaxa(Fpxdel,az,4);sf_oaxa(Fpxdel,ax,5);sf_oaxa(Fpxdel,ay,6);
	sf_oaxa(Fpydel,az,4);sf_oaxa(Fpydel,ax,5);sf_oaxa(Fpydel,ay,6);

	sf_oaxa(Fvzdel,az,4);sf_oaxa(Fvzdel,ax,5);sf_oaxa(Fvzdel,ay,6);
	sf_oaxa(Fvxdel,az,4);sf_oaxa(Fvxdel,ax,5);sf_oaxa(Fvxdel,ay,6);
	sf_oaxa(Fvydel,az,4);sf_oaxa(Fvydel,ax,5);sf_oaxa(Fvydel,ay,6);

	sf_oaxa(Fhzdel,az,4);sf_oaxa(Fhzdel,ax,5);sf_oaxa(Fhzdel,ay,6);
	sf_oaxa(Fhxdel,az,4);sf_oaxa(Fhxdel,ax,5);sf_oaxa(Fhxdel,ay,6);
	sf_oaxa(Fhydel,az,4);sf_oaxa(Fhydel,ax,5);sf_oaxa(Fhydel,ay,6);



    }

    /*------------------------------------------------------------*/
    
    /*------------------------------------------------------------*/
    /* read model params */
    c11=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c12=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c13=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c22=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c23=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c33=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c44=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c55=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c66=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    /*   tmp=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay)); */
    
    nnn = sf_n(ax)*sf_n(ay)*sf_n(az);
    /*  sf_floatread(c11[0][0],nnn,Fccc);  */
/*     sf_floatread(c12[0][0],nnn,Fccc); */
/*     sf_floatread(c13[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(c22[0][0],nnn,Fccc); */
/*     sf_floatread(c23[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(c33[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(c44[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(c55[0][0],nnn,Fccc); */
/*     sf_floatread(tmp[0][0],nnn,Fccc); */
/*     sf_floatread(c66[0][0],nnn,Fccc); */

    sf_floatread(c11[0][0],nnn,Fccc); 
    sf_floatread(c22[0][0],nnn,Fccc);
    sf_floatread(c33[0][0],nnn,Fccc);
    sf_floatread(c44[0][0],nnn,Fccc);
    sf_floatread(c55[0][0],nnn,Fccc);
    sf_floatread(c66[0][0],nnn,Fccc);
    sf_floatread(c12[0][0],nnn,Fccc);
    sf_floatread(c13[0][0],nnn,Fccc);
    sf_floatread(c23[0][0],nnn,Fccc);
    /*------------------------------------------------------------*/
    /* compute derivatives */
    
   
    
    if(stat) {
	if(domain[0]=='k'){ 
	    sf_warning("stop0");
	    wfsepK3(pxdel,pydel,pzdel,
		    vxdel,vydel,vzdel,
		    hxdel,hydel,hzdel,
		    ax,ay,az,
		    c11[0][0][0],c12[0][0][0],c13[0][0][0],
		    c22[0][0][0],c23[0][0][0],c33[0][0][0],
		    c44[0][0][0],c55[0][0][0],c66[0][0][0],
		    wfs3);
	    ndel=sf_n(ax)*sf_n(ay)*sf_n(az);
	    sf_floatwrite(pzdel[0][0],ndel,Fpzdel);
	    sf_floatwrite(pxdel[0][0],ndel,Fpxdel);
	    sf_floatwrite(pydel[0][0],ndel,Fpydel);

	    sf_floatwrite(vzdel[0][0],ndel,Fvzdel);
	    sf_floatwrite(vxdel[0][0],ndel,Fvxdel);
	    sf_floatwrite(vydel[0][0],ndel,Fvydel);

	    sf_floatwrite(hzdel[0][0],ndel,Fhzdel);
	    sf_floatwrite(hxdel[0][0],ndel,Fhxdel);
	    sf_floatwrite(hydel[0][0],ndel,Fhydel);
	    
	}
	else{
	    sf_warning("stop1");
	    wfsep3(pxdel,pydel,pzdel,
		   vxdel,vydel,vzdel,
		   hxdel,hydel,hzdel,
		   sx,sy,sz,
		   c11[0][0][0],c12[0][0][0],c13[0][0][0],
		   c22[0][0][0],c23[0][0][0],c33[0][0][0],
		   c44[0][0][0],c55[0][0][0],c66[0][0][0],
		   wfs3);
	    ndel=sf_n(sx)*sf_n(sy)*sf_n(sz);
	    sf_floatwrite(pzdel[0][0],ndel,Fpzdel);
	    sf_floatwrite(pxdel[0][0],ndel,Fpxdel);
	    sf_floatwrite(pydel[0][0],ndel,Fpydel);
	    
	    sf_floatwrite(vzdel[0][0],ndel,Fvzdel);
	    sf_floatwrite(vxdel[0][0],ndel,Fvxdel);
	    sf_floatwrite(vydel[0][0],ndel,Fvydel);

	    sf_floatwrite(hzdel[0][0],ndel,Fhzdel);
	    sf_floatwrite(hxdel[0][0],ndel,Fhxdel);
	    sf_floatwrite(hydel[0][0],ndel,Fhydel);
	    
	}

    }
    
    
		
    else {

	
	if(domain[0]=='k'){
	    for(    iy=0;iy<sf_n(ay);iy++){
		for(    ix=0;ix<sf_n(ax);ix++){
		    for(iz=0;iz<sf_n(az);iz++){
			if(verb) fprintf(stderr,"%5d %5d %5d",iz,ix, iy);
			wfsepK3(pxdel,pydel,pzdel,
				vxdel,vydel,vzdel,
				hxdel,hydel,hzdel,
				ax,ay,az,
				c11[iy][ix][iz],c12[iy][ix][iz],c13[iy][ix][iz],
				c22[iy][ix][iz],c23[iy][ix][iz],c33[iy][ix][iz],
				c44[iy][ix][iz],c55[iy][ix][iz],c66[iy][ix][iz],
				wfs3);
			ndel=sf_n(ax)*sf_n(ay)*sf_n(az);
			sf_floatwrite(pzdel[0][0],ndel,Fpzdel);
			sf_floatwrite(pxdel[0][0],ndel,Fpxdel);
			sf_floatwrite(pydel[0][0],ndel,Fpydel);
			
			sf_floatwrite(vzdel[0][0],ndel,Fvzdel);
			sf_floatwrite(vxdel[0][0],ndel,Fvxdel);
			sf_floatwrite(vydel[0][0],ndel,Fvydel);
			
			sf_floatwrite(hzdel[0][0],ndel,Fhzdel);
			sf_floatwrite(hxdel[0][0],ndel,Fhxdel);
			sf_floatwrite(hydel[0][0],ndel,Fhydel);
			
			if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b ");
		    }   
		}
	    }
	}
	else{    
	    for(    iy=0;iy<sf_n(ay);iy++){
		for(    ix=0;ix<sf_n(ax);ix++){
		    for(iz=0;iz<sf_n(az);iz++){
			
			if(verb) fprintf(stderr,"%5d %5d %5d",iz,ix,iy);
			wfsep3(pxdel,pydel,pzdel,
			       vxdel,vydel,vzdel,
			       hxdel,hydel,hzdel,
			       sx,sy,sz,
			       c11[iy][ix][iz],c12[iy][ix][iz],c13[iy][ix][iz],
			       c22[iy][ix][iz],c23[iy][ix][iz],c33[iy][ix][iz],
			       c44[iy][ix][iz],c55[iy][ix][iz],c66[iy][ix][iz],
			       wfs3);
			ndel=sf_n(sx)*sf_n(sy)*sf_n(sz);
			sf_floatwrite(pzdel[0][0],ndel,Fpzdel);
			sf_floatwrite(pxdel[0][0],ndel,Fpxdel);
			sf_floatwrite(pydel[0][0],ndel,Fpydel);


			sf_floatwrite(vzdel[0][0],ndel,Fvzdel);
			sf_floatwrite(vxdel[0][0],ndel,Fvxdel);
			sf_floatwrite(vydel[0][0],ndel,Fvydel);
			
			sf_floatwrite(hzdel[0][0],ndel,Fhzdel);
			sf_floatwrite(hxdel[0][0],ndel,Fhxdel);
			sf_floatwrite(hydel[0][0],ndel,Fhydel);
			



			if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b ");
		    }
		}	    
	    }
	}
	/*------------------------------------------------------------*/

	wfsep3_close(wfs3);
	exit (0);
    }
    
    
}
