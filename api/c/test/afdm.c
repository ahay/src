/* time-domain acoustic FD modeling */
#include <rsf.h>
int main(int argc, char* argv[])
{
    /* Laplacian coefficients */
    float c0=-30./12.,c1=+16./12.,c2=- 1./12.;

    bool verb;           /* verbose flag */
    sf_file Fw=NULL,Fv=NULL,Fr=NULL,Fo=NULL; /* I/O files */
    sf_axis at,az,ax;    /* cube axes */
    int it,iz,ix;        /* index variables */
    int nt,nz,nx;
    float dt,dz,dx,idx,idz,dt2;

    float  *ww,**vv,**rr;     /* I/O arrays*/
    float **um,**uo,**up,**ud;/* tmp arrays */

    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=0; /* verbose flag */

    /* setup I/O files */
    Fw = sf_input ("in" );
    Fo = sf_output("out");
    Fv = sf_input ("vel");
    Fr = sf_input ("ref");

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fv,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); dx = sf_d(ax);

    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);

    dt2 =    dt*dt;
    idz = 1/(dz*dz);
    idx = 1/(dx*dx);

    /* read wavelet, velocity & reflectivity */
    ww=sf_floatalloc(nt);     sf_floatread(ww   ,nt   ,Fw);
    vv=sf_floatalloc2(nz,nx); sf_floatread(vv[0],nz*nx,Fv);
    rr=sf_floatalloc2(nz,nx); sf_floatread(rr[0],nz*nx,Fr);

    /* allocate temporary arrays */
    um=sf_floatalloc2(nz,nx);
    uo=sf_floatalloc2(nz,nx);
    up=sf_floatalloc2(nz,nx);
    ud=sf_floatalloc2(nz,nx);

    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ud[ix][iz]=0;
	}
    }

    /* MAIN LOOP */
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	/* 4th order laplacian */
	for (iz=2; iz<nz-2; iz++) {
	    for (ix=2; ix<nx-2; ix++) {
		ud[ix][iz] = 
		    c0* uo[ix  ][iz  ] * (idx+idz) + 
		    c1*(uo[ix-1][iz  ] + uo[ix+1][iz  ])*idx +
		    c2*(uo[ix-2][iz  ] + uo[ix+2][iz  ])*idx +
		    c1*(uo[ix  ][iz-1] + uo[ix  ][iz+1])*idz +
		    c2*(uo[ix  ][iz-2] + uo[ix  ][iz+2])*idz;	  
	    }
	}

	/* inject wavelet */
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ud[ix][iz] -= ww[it] * rr[ix][iz];
	    }
	}

	/* scale by velocity */
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ud[ix][iz] *= vv[ix][iz]*vv[ix][iz];
	    }
	}
	
	/* time step */
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		up[ix][iz] = 
		    2*uo[ix][iz] 
		    - um[ix][iz] 
		    + ud[ix][iz] * dt2; 
		
		um[ix][iz] = uo[ix][iz];
		uo[ix][iz] = up[ix][iz];
	    }
	}
	
	/* write wavefield to output */
	sf_floatwrite(uo[0],nz*nx,Fo);
    }
    if(verb) fprintf(stderr,"\n");

    exit (0);
}
