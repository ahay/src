/* 
 * time-domain acoustic FD modeling
 */

#include <rsf.h>
/*#include "cub.h"*/

int main(int argc, char* argv[])
{
    bool verb; /* verbose flag */

    /* I/O files */
    sf_file Fw,Fv,Fr;
    sf_file Fo;

    /* cube axes */
    axa at,az,ax;
    int it,iz,ix;
    float idx,idz,dt2;

    /* arrays */
    float  *ww, **vv, **rr;    
    float **um, **uo, **up, **ud;

    /* Laplacian coefficients */
    float c0=-30./12.;
    float c1=+16./12.;
    float c2=- 1./12.;

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=0;

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fv = sf_input ("vel");
    Fr = sf_input ("ref");

    /* read axes*/
    iaxa(Fw,&at,1);
    iaxa(Fv,&az,1);
    iaxa(Fv,&ax,2);

    /* setup output header */
    oaxa(Fo,&az,1);
    oaxa(Fo,&ax,2);
    oaxa(Fo,&at,3);

    dt2 =    at.d*at.d;
    idz = 1/(az.d*az.d);
    idx = 1/(ax.d*ax.d);
     
    /* read wavelet, velocity & reflectivity */
    ww=sf_floatalloc(at.n);
    vv=sf_floatalloc2(az.n,ax.n);
    rr=sf_floatalloc2(az.n,ax.n);

    sf_floatread(ww   ,at.n     ,Fw);    
    sf_floatread(vv[0],az.n*ax.n,Fv);
    sf_floatread(rr[0],az.n*ax.n,Fr);

    /* allocate temporary arrays */
    um=sf_floatalloc2(az.n,ax.n);
    uo=sf_floatalloc2(az.n,ax.n);
    up=sf_floatalloc2(az.n,ax.n);
    ud=sf_floatalloc2(az.n,ax.n);
    
    for (iz=0; iz<az.n; iz++) {
	for (ix=0; ix<ax.n; ix++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ud[ix][iz]=0;
	}
    }
    
    /* 
     *  MAIN LOOP
     */
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<at.n; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	/* 4th order laplacian */
	for (iz=2; iz<az.n-2; iz++) {
	    for (ix=2; ix<ax.n-2; ix++) {
		ud[ix][iz] = 
		    c0* uo[ix  ][iz  ] * (idx+idz) + 
		    c1*(uo[ix-1][iz  ] + uo[ix+1][iz  ])*idx +
		    c2*(uo[ix-2][iz  ] + uo[ix+2][iz  ])*idx +
		    c1*(uo[ix  ][iz-1] + uo[ix  ][iz+1])*idz +
		    c2*(uo[ix  ][iz-2] + uo[ix  ][iz+2])*idz;	  
	    }
	}
	
	for (iz=0; iz<az.n; iz++) {
	    for (ix=0; ix<ax.n; ix++) {
		ud[ix][iz] *= vv[ix][iz]*vv[ix][iz];
	    }
	}
	
	/* inject wavelet */
	for (iz=0; iz<az.n; iz++) {
	    for (ix=0; ix<ax.n; ix++) {
		ud[ix][iz] -= ww[it] * rr[ix][iz];
	    }
	}

	/* time step */
	for (iz=0; iz<az.n; iz++) {
	    for (ix=0; ix<ax.n; ix++) {
		up[ix][iz] = 
		    2*uo[ix][iz] 
		    - um[ix][iz] 
		    + ud[ix][iz] * dt2; 
		
		um[ix][iz] = uo[ix][iz];
		uo[ix][iz] = up[ix][iz];
	    }
	}
	
	/* write wavefield to output */
	sf_floatwrite(uo[0],az.n*ax.n,Fo);
    }
    if(verb) fprintf(stderr,"\n");    

    exit (0);
}
