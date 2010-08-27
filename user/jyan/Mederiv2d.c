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
#include "tool2D.h"

/*^*/

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/*------------------------------------------------------------*/

#define  KMAP(i,n) (i<n/2.) ? SF_PI*i/(n/2.) : SF_PI*(-n+i)/(n/2.);
#define  KMAPK(i,n) (  ( i-n/2 )/(n/2.)*SF_PI  );
/*#define  KMAP(i,n) (  ( i-n/2 )/(n/2.)*SF_PI  );*/

/* Joe's taper */
/*#define TAPER(k) (0.5*(1+cosf(k))) */

/* 2nd order */
#define TAPER2(k) (k!=0 ? sin(  k)/k : 1)
/*#define TAPER(k) (sin(  k)/k )*/
/* 4th order */
#define TAPER4(k) (k!=0 ?	       \
		  4./3.  *sin(  k)/k - \
		  1./6.  *sin(2*k)/k : 1)

/* 6th order */
#define TAPER6(k) (k!=0 ?	       \
		  3./2.  *sin(  k)/k - \
		  3./10. *sin(2*k)/k + \
		  1./60. *sin(3*k)/k : 1)

/* 8th order */
#define TAPER8(k) (k!=0 ?			\
		  8./5.  *sin(  k)/(k) -	\
		  2./5.  *sin(2*k)/(k) +	\
		  8./105.*sin(3*k)/(k) -	\
		  1./140.*sin(4*k)/(k) : 1)

/*#define TAPER(k) 1.*/

/* Dave's 2nd order */
#define TAPERD(k1,k2) (k1!=0 ? sinf(k1/2)*cosf(k2/2)/k1 :  2*cosf(k2/2) )

#define TAPERG(kmag,sig) exp(-(kmag*kmag)/2/sig/sig)/2.


#define TAPERS(k) (k!=0 ? sin(  k)/k : 1)

/*#define filter(k) sin(k)/k*/
#define TAPER(k) (k!=0 ?			\
		  (k<SF_PI ? filter(k) :0)	\
		  :1)

#define filter(k)	  8./5.  *sin(  k)/k -			\
					  2./5.  *sin(2*k)/k +	\
					  8./105.*sin(3*k)/k -  \
					  1./140.*sin(4*k)/k 
#define TAPERK(k1,k2) (k1*k1+k2*k2 != 0 ?				\
		       (sqrt(k1*k1+k2*k2)<SF_PI  ? filter(k1) :0.)	\
		       :1)

/*Butterworth 3rd order low pass filter*/
#define filterB(k)    (1+cos(k))*(1+cos(k))*(1+cos(k))/(5+3*cos(2*k))
#define filterBreal(k)   4*pow(cos(k/2),4) * (-1+3*cos(k))         /(5+3*cos(2*k))
#define filterBimag(k)   4*pow(cos(k/2),3) * ( 1+3*cos(k))*sin(k/2)/(5+3*cos(2*k))

#define TAPERB(kmag)  (kmag<SF_PI  ? filterB(kmag)  :0 )

#define TAPERBR(kmag) (kmag<SF_PI  ? filterBreal(kmag)  :0 )

#define TAPERBI(kmag) (kmag<SF_PI  ? filterBimag(kmag)  :0 ) 
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/

float **tpvalue;
char *domain;
char *tapertype;

int main(int argc, char* argv[])
{
    bool verb;
    bool stat;

  
    
    /* I/O files */
    sf_file Fccc=NULL;

    sf_file Fspk=NULL;  /* template for derivative operator */    
    sf_file Fzdel=NULL; /* z derivative */
    sf_file Fxdel=NULL; /* x derivative */
    
    sf_axis ax,az;      /* cube axes */
    sf_axis sx,sz;      /* cube axes */
    
    float **c11,**c13,**c15,**c33,**c35,**c55; /* input file */
    float  mc11, mc13, mc15, mc33, mc35, mc55;
    float **xdel,**zdel;                /* output derivative */
    
    int ix,iz,nn;
    wfs2d wfs;
    int order;
    float sig;
/*    float sig;*/

    float k2;
    int nx,nz;
    int jx,jz;
    float kx,kz;
    
    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("stat",&stat)) stat=true;  /* stationary operator */
    
/*    sf_warning("\n stat= %d \n",stat);*/
    
    if(NULL== (domain=sf_getstring("domain"))) domain="k";    
    if(NULL== (tapertype=sf_getstring("tapertype"))) tapertype="g";
    if(! sf_getfloat("sig",&sig)) sig=1.5;  /* sigma */    
    if(! sf_getint("order",&order)) order=8;  /* order */


    if(verb) sf_warning("domain= %s",domain);
    if(verb) sf_warning("tapertype= %s",tapertype);
    if(verb && tapertype[0]=='s') sf_warning("order= %d",order);
    if(verb && tapertype[0]=='g') sf_warning("sig= %f",sig);
    /*------------------------------------------------------------*/
    /* input files and axes */
    Fccc = sf_input ("ccc");      
    az = sf_iaxa(Fccc,1);
    ax = sf_iaxa(Fccc,2);

    Fspk = sf_input ("in");
    sz = sf_iaxa(Fspk,1);
    sx = sf_iaxa(Fspk,2);
  
    if(verb){
	sf_raxa(az);
	sf_raxa(ax);
	sf_raxa(sz);
	sf_raxa(sx);
	sf_warning("");
    }

    /*------------------------------------------------------------*/
    /* output files*/
    Fzdel = sf_output ("zdel");
    Fxdel = sf_output ("xdel");
    if(domain[0]=='k'){
	sf_oaxa(Fxdel,sz,1);
	sf_oaxa(Fxdel,sx,2);
	sf_oaxa(Fzdel,sz,1);
	sf_oaxa(Fzdel,sx,2);
	zdel=sf_floatalloc2(sf_n(sz),sf_n(sx));
	xdel=sf_floatalloc2(sf_n(sz),sf_n(sx));
	wfs = wfsep_init(sx,sz);
    }
    else{
	sf_oaxa(Fxdel,sz,1);
	sf_oaxa(Fxdel,sx,2);
	sf_oaxa(Fzdel,sz,1);
	sf_oaxa(Fzdel,sx,2);
	zdel=sf_floatalloc2(sf_n(sz),sf_n(sx));
	xdel=sf_floatalloc2(sf_n(sz),sf_n(sx));
	wfs = wfsep_init(sx,sz);
    }
        
    if(!stat) {
	sf_oaxa(Fzdel,az,3);
	sf_oaxa(Fzdel,ax,4);
	sf_oaxa(Fxdel,az,3);
	sf_oaxa(Fxdel,ax,4);
    }

    

    /*------------------------------------------------------------*/
    /* read model params */
    c11=sf_floatalloc2(sf_n(az),sf_n(ax));
    c13=sf_floatalloc2(sf_n(az),sf_n(ax));
    c15=sf_floatalloc2(sf_n(az),sf_n(ax));
    c33=sf_floatalloc2(sf_n(az),sf_n(ax));
    c35=sf_floatalloc2(sf_n(az),sf_n(ax));
    c55=sf_floatalloc2(sf_n(az),sf_n(ax));

    
    nn = sf_n(az)*sf_n(ax);
    sf_floatread(c11[0],nn,Fccc); 
    sf_floatread(c13[0],nn,Fccc);
    sf_floatread(c15[0],nn,Fccc);
    sf_floatread(c33[0],nn,Fccc);
    sf_floatread(c35[0],nn,Fccc);
    sf_floatread(c55[0],nn,Fccc);



    
    /*------------------------------------------------------------*/
    /* allocate arrays for derivatives */
    if (domain[0]=='k'){
	nx=sf_n(sx);
	nz=sf_n(sz);
	sf_warning("nx=%d nz=%d",nx,nz);
        tpvalue=sf_floatalloc2(nz,nx);
	if(tapertype[0]=='g'){ 
	    for(    jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
		for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		    k2=sqrt(kx*kx+kz*kz);
		    tpvalue[jx][jz]=TAPERG(k2,sig);
		}
	    }
	}
	else if(tapertype[0]=='s'){
	    for(    jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
		for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		    switch(order){
			case 8: tpvalue[jx][jz]=TAPER8(kx); break;
			case 6: tpvalue[jx][jz]=TAPER6(kx); break;
			case 4: tpvalue[jx][jz]=TAPER4(kx); break;
			case 2: tpvalue[jx][jz]=TAPER2(kx); 
		    }
		}
	    }
	}
	else{
	    for(    jx=0;jx<nx;jx++){
		for(jz=0;jz<nz;jz++){
		    tpvalue[jx][jz]=1.0;
		}
	    }
	}
    }
    else{
	nx=sf_n(sx);
	nz=sf_n(sz);
	sf_warning("nx=%d nz=%d",nx,nz);tpvalue=sf_floatalloc2(nz,nx);
	if(tapertype[0]=='g'){ 
	    for(    jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
		for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		    k2=sqrt(kx*kx+kz*kz);
		    tpvalue[jx][jz]=TAPERG(k2,sig);
		}
	    }
	}
	else if(tapertype[0]=='s'){
	    for(    jx=0;jx<nx;jx++){kx = KMAP(jx,nx);
		for(jz=0;jz<nz;jz++){kz = KMAP(jz,nz);
		    switch(order){
			case 8: tpvalue[jx][jz]=TAPER8(kx); break;
			case 6: tpvalue[jx][jz]=TAPER6(kx); break;
			case 4: tpvalue[jx][jz]=TAPER4(kx); break;
			case 2: tpvalue[jx][jz]=TAPER2(kx); 
		    }
		}
	    }
	}
	else{
	    for(    jx=0;jx<nx;jx++){
		for(jz=0;jz<nz;jz++){
		    tpvalue[jx][jz]=1.0;
		}
	    }
	}
    }




    /*------------------------------------------------------------*/
    /* compute derivatives */
    
   
    
    if(stat) {
	mc11 = sf_quantile(nn/2,nn,c11[0]);
	mc13 = sf_quantile(nn/2,nn,c13[0]);
	mc15 = sf_quantile(nn/2,nn,c15[0]);
	mc33 = sf_quantile(nn/2,nn,c33[0]);
	mc35 = sf_quantile(nn/2,nn,c35[0]);
	mc55 = sf_quantile(nn/2,nn,c55[0]);
	
	if(verb) {
	    sf_warning("c11=%g c13=%g c15=%g c33=%g c35=%g c55=%g",
		       mc11,  mc13,  mc15,  mc33,  mc35,  mc55);
	}
	
	
	if(domain[0]=='k'){ 
	    wfsep(zdel,xdel,sx,sz,
		   mc11,mc13,mc15,mc33,mc35,mc55,
		   wfs);
	    sf_floatwrite(zdel[0],sf_n(sz)*sf_n(sx),Fzdel);
	    sf_floatwrite(xdel[0],sf_n(sz)*sf_n(sx),Fxdel);
	}
	else{
	    wfsep(zdel,xdel,sx,sz,
		  mc11,mc13,mc15,mc33,mc35,mc55,
		  wfs);
	    sf_floatwrite(zdel[0],sf_n(sz)*sf_n(sx),Fzdel);
	    sf_floatwrite(xdel[0],sf_n(sz)*sf_n(sx),Fxdel);
	}

    }
    

		
    else {

	
	if(domain[0]=='k'){
	    for(    ix=0;ix<sf_n(ax);ix++){
		for(iz=0;iz<sf_n(az);iz++){
		    if(verb) fprintf(stderr,"%5d %5d",iz,ix);
		    
		    wfsep(zdel,xdel,ax,az,
			   c11[ix][iz],c13[ix][iz],c15[ix][iz],c33[ix][iz],
			   c35[ix][iz],c55[ix][iz],wfs);
		    sf_floatwrite(zdel[0],sf_n(az)*sf_n(ax),Fzdel);
		    sf_floatwrite(xdel[0],sf_n(az)*sf_n(ax),Fxdel);
		    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b ");
		}   
	    }
	}
	else{    
	    
	    for(    ix=0;ix<sf_n(ax);ix++){
		for(iz=0;iz<sf_n(az);iz++){
		    
		    if(verb) fprintf(stderr,"%5d %5d",iz,ix);
		    wfsep(zdel,xdel,sx,sz,
			  c11[ix][iz],c13[ix][iz],c15[ix][iz],c33[ix][iz],
			  c35[ix][iz],c55[ix][iz],wfs);
		    
		    sf_floatwrite(zdel[0],sf_n(sz)*sf_n(sx),Fzdel);
		    sf_floatwrite(xdel[0],sf_n(sz)*sf_n(sx),Fxdel);
		    
		    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b ");
		}
	    }	    
	}
	/*------------------------------------------------------------*/
	wfsep_close(wfs);
	
	exit (0);
    }
    
    
}
