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
#include "tool3tti.h"
#include "dsyevv3.h"
#include "dsyevc3.h"

/* These are the programs to calculte 3x3 eigenvalues and eigenvectors from Kopp */
/* 's 2008 paper:  */
/* Efficient numerical diagonalization of hermitian 3x3
 * matrices */
/* International journal of modern physics C Vol 19, No 3 523-548 */
/* downloaded from www.mpi-hd.mpg.de/~globes/3x3  */
#include <math.h>


int main(int argc, char* argv[])
{
    bool verb;
    bool stat;

    char *domain;
    
    /* I/O files */
    sf_file Fccc=NULL;      
       
    sf_file Fspk=NULL;  /* template for derivative operator */    
    sf_file Fzdel=NULL; /* z derivative */
    sf_file Fxdel=NULL; /* x derivative */
    sf_file Fydel=NULL; /* y derivative */


    sf_axis ax,ay,az;      /* cube axes */
    sf_axis sx,sy,sz;      /* cube axes */
    
    
    float ***c11,***c12,***c13,***c14,***c15,***c16;
    float ***c22,***c23,***c24,***c25,***c26;
    float ***c33,***c34,***c35,***c36;
    float ***c44,***c45,***c46;
    float ***c55,***c56;
    float ***c66;
    

    float ***xdel,***ydel,***zdel;                /* output derivative */
    
    int ix,iy,iz,nnn,ndel;
    wfs3d wfs3;
    
    /*------------------------------------------------------------*/
   



    /*------------------------------------------------------------*/

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
    Fzdel = sf_output ("zdel");
    Fxdel = sf_output ("xdel");
    Fydel = sf_output ("ydel");
    if(domain[0]=='k'){
	sf_oaxa(Fxdel,az,1);sf_oaxa(Fxdel,ax,2);sf_oaxa(Fxdel,ay,3);
	sf_oaxa(Fydel,az,1);sf_oaxa(Fydel,ax,2);sf_oaxa(Fydel,ay,3);	
	sf_oaxa(Fzdel,az,1);sf_oaxa(Fzdel,ax,2);sf_oaxa(Fzdel,ay,3);
	zdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	ydel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	xdel=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
	wfs3 = wfsep3_init(ax,ay,az);
    }
    else{
	sf_oaxa(Fxdel,sz,1);sf_oaxa(Fxdel,sx,2);sf_oaxa(Fxdel,sy,3);
	sf_oaxa(Fydel,sz,1);sf_oaxa(Fydel,sx,2);sf_oaxa(Fydel,sy,3);
	sf_oaxa(Fzdel,sz,1);sf_oaxa(Fzdel,sx,2);sf_oaxa(Fzdel,sy,3);
	zdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	ydel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	xdel=sf_floatalloc3(sf_n(sz),sf_n(sx),sf_n(sy));
	wfs3 = wfsep3_init(sx,sy,sz);
    }
        
    if(!stat) {
	sf_oaxa(Fzdel,az,4);sf_oaxa(Fzdel,ax,5);sf_oaxa(Fzdel,ay,6);
	sf_oaxa(Fxdel,az,4);sf_oaxa(Fxdel,ax,5);sf_oaxa(Fxdel,ay,6);
	sf_oaxa(Fydel,az,4);sf_oaxa(Fydel,ax,5);sf_oaxa(Fydel,ay,6);
    }

    /*------------------------------------------------------------*/
    
    /*------------------------------------------------------------*/
    /* read model params */
    c11=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c12=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c13=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c14=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c15=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c16=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c22=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c23=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c24=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c25=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c26=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c33=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c34=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c35=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c36=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c44=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c45=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c46=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c55=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c56=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
    c66=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));
   
    

    nnn = sf_n(ax)*sf_n(ay)*sf_n(az);
    sf_floatread(c11[0][0],nnn,Fccc);
    sf_floatread(c12[0][0],nnn,Fccc);
    sf_floatread(c13[0][0],nnn,Fccc);
    sf_floatread(c14[0][0],nnn,Fccc);
    sf_floatread(c15[0][0],nnn,Fccc);
    sf_floatread(c16[0][0],nnn,Fccc);
    sf_floatread(c22[0][0],nnn,Fccc);
    sf_floatread(c23[0][0],nnn,Fccc);
    sf_floatread(c24[0][0],nnn,Fccc);
    sf_floatread(c25[0][0],nnn,Fccc);
    sf_floatread(c26[0][0],nnn,Fccc);
    sf_floatread(c33[0][0],nnn,Fccc);
    sf_floatread(c34[0][0],nnn,Fccc);
    sf_floatread(c35[0][0],nnn,Fccc);
    sf_floatread(c36[0][0],nnn,Fccc);
    sf_floatread(c44[0][0],nnn,Fccc);
    sf_floatread(c45[0][0],nnn,Fccc);
    sf_floatread(c46[0][0],nnn,Fccc);
    sf_floatread(c55[0][0],nnn,Fccc);
    sf_floatread(c56[0][0],nnn,Fccc);
    sf_floatread(c66[0][0],nnn,Fccc);
  
    /*------------------------------------------------------------*/
    /* compute derivatives */
    sf_warning("stiffness matrix");
    sf_warning("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", 
	       c11[0][0][0],c12[0][0][0],c13[0][0][0],c14[0][0][0],c15[0][0][0],c16[0][0][0]);
    sf_warning("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", 
 	       c12[0][0][0],c22[0][0][0],c23[0][0][0],c24[0][0][0],c25[0][0][0],c26[0][0][0]);
    sf_warning("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", 
	       c13[0][0][0],c23[0][0][0],c33[0][0][0],c34[0][0][0],c35[0][0][0],c36[0][0][0]);
    sf_warning("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", 
	       c14[0][0][0],c24[0][0][0],c34[0][0][0],c44[0][0][0],c45[0][0][0],c46[0][0][0]);
    sf_warning("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", 
	       c15[0][0][0],c25[0][0][0],c35[0][0][0],c45[0][0][0],c55[0][0][0],c56[0][0][0]);
    sf_warning("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", 
	       c16[0][0][0],c26[0][0][0],c36[0][0][0],c46[0][0][0],c56[0][0][0],c66[0][0][0]);
	
    if(stat) {
		
	if(domain[0]=='k'){ 
	    sf_warning("stop0");
	    wfsepK3(xdel,ydel,zdel,ax,ay,az,		        
		    c11[0][0][0],c12[0][0][0],c13[0][0][0],c14[0][0][0],c15[0][0][0],c16[0][0][0],
		    c22[0][0][0],c23[0][0][0],c24[0][0][0],c25[0][0][0],c26[0][0][0],
		    c33[0][0][0],c34[0][0][0],c35[0][0][0],c36[0][0][0],
		    c44[0][0][0],c45[0][0][0],c46[0][0][0],
		    c55[0][0][0],c56[0][0][0],
		    c66[0][0][0],
		    wfs3);
	    ndel=sf_n(ax)*sf_n(ay)*sf_n(az);
	    sf_floatwrite(zdel[0][0],ndel,Fzdel);
	    sf_floatwrite(xdel[0][0],ndel,Fxdel);
	    sf_floatwrite(ydel[0][0],ndel,Fydel);
	}
	else{
	    wfsep3(xdel,ydel,zdel,sx,sy,sz,		        
		    c11[0][0][0],c12[0][0][0],c13[0][0][0],c14[0][0][0],c15[0][0][0],c16[0][0][0],
		    c22[0][0][0],c23[0][0][0],c24[0][0][0],c25[0][0][0],c26[0][0][0],
		    c33[0][0][0],c34[0][0][0],c35[0][0][0],c36[0][0][0],
		    c44[0][0][0],c45[0][0][0],c46[0][0][0],
		    c55[0][0][0],c56[0][0][0],
		    c66[0][0][0],
		    wfs3);
	    ndel=sf_n(sx)*sf_n(sy)*sf_n(sz);
	    sf_floatwrite(zdel[0][0],ndel,Fzdel);
	    sf_floatwrite(xdel[0][0],ndel,Fxdel);
	    sf_floatwrite(ydel[0][0],ndel,Fydel);
	}

    }
    
    
		
    else {

	
	if(domain[0]=='k'){
	    for(    iy=0;iy<sf_n(ay);iy++){
		for(    ix=0;ix<sf_n(ax);ix++){
		    for(iz=0;iz<sf_n(az);iz++){
			if(verb) fprintf(stderr,"%5d %5d %5d",iz,ix, iy);
			wfsepK3(xdel,ydel,zdel,ax,ay,az,		        
				c11[0][0][0],c12[0][0][0],c13[0][0][0],c14[0][0][0],c15[0][0][0],c16[0][0][0],
				c22[0][0][0],c23[0][0][0],c24[0][0][0],c25[0][0][0],c26[0][0][0],
				c33[0][0][0],c34[0][0][0],c35[0][0][0],c36[0][0][0],
				c44[0][0][0],c45[0][0][0],c46[0][0][0],
				c55[0][0][0],c56[0][0][0],
				c66[0][0][0],
				wfs3);
			
			ndel=sf_n(ax)*sf_n(ay)*sf_n(az);
			sf_floatwrite(zdel[0][0],ndel,Fzdel);
			sf_floatwrite(xdel[0][0],ndel,Fxdel);
			sf_floatwrite(ydel[0][0],ndel,Fydel);
			if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b ");
		    }   
		}
	    }
	}
	else{    
	    for(    iy=0;iy<sf_n(ay);iy++){
		for(    ix=0;ix<sf_n(ax);ix++){
		    for(iz=0;iz<sf_n(az);iz++){
			
			if(verb) fprintf(stderr,"%5d %5d %5d",iz,ix, iy);
			wfsep3(xdel,ydel,zdel,ax,ay,az,		        
				c11[0][0][0],c12[0][0][0],c13[0][0][0],c14[0][0][0],c15[0][0][0],c16[0][0][0],
				c22[0][0][0],c23[0][0][0],c24[0][0][0],c25[0][0][0],c26[0][0][0],
				c33[0][0][0],c34[0][0][0],c35[0][0][0],c36[0][0][0],
				c44[0][0][0],c45[0][0][0],c46[0][0][0],
				c55[0][0][0],c56[0][0][0],
				c66[0][0][0],
				wfs3);
			ndel=sf_n(sx)*sf_n(sy)*sf_n(sz);
			sf_floatwrite(zdel[0][0],ndel,Fzdel);
			sf_floatwrite(xdel[0][0],ndel,Fxdel);
			sf_floatwrite(ydel[0][0],ndel,Fydel);
			if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b ");
		    }
		}	    
	    }
	}
	/*------------------------------------------------------------*/
	sf_warning("stop1");
	wfsep3_close(wfs3);
	sf_warning("stop1");
	exit (0);
    }
    
    
}
