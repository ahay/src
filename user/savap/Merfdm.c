
/* 
 * time-domain acoustic FD modeling
 */


#include <rsf.h>

int main(int argc, char* argv[])
{

    /* I/O files */
    sf_file in,out;
    sf_file vel; /* velocity     */
    sf_file den; /* density      */
    sf_file ref; /* reflectivity */

    /* I/O arrays */
    float  *wvl;
    float **vv;
    float **dd;
    float **rr;
    
    /* cube axes */
    int   nt,nz,nx;
    float ot,oz,ox;
    float dt,dz,dx;
    float idx,idz;

    /* wavefield arrays */
    float **um, **uo, **up, **ud;

    /* Laplacian coefficients */
    float c0=-30./12.;
    float c1=+16./12.;
    float c2=- 1./12.;

    /* misc */
    int it,iz,ix;

    
    /* 
     * init param I/O
     */
    sf_init(argc,argv);

    /* 
     * setup I/O files
     */
    in  = sf_input ("in" );
    out = sf_output("out");
    vel = sf_input ("vel");
    den = sf_input ("den");
    ref = sf_input ("ref");
     
    /* 
     * read wavelet params
     */
    if(! sf_histint(   in,"n1",&nt)) sf_error("No nt in wvl");  
    if(! sf_histfloat( in,"d1",&dt)) sf_error("No dt in wvl");
    if(! sf_histfloat( in,"o1",&ot)) sf_error("No ot in wvl");
    
    /* 
     * read wavelet
     */
    wvl=sf_floatalloc(nt);
    sf_floatread(wvl,nt,in);    
    
    /* 
     * read model params
     */
    if(! sf_histint(  vel,"n1",&nz)) sf_error("No nz in vel");  
    if(! sf_histfloat(vel,"d1",&dz)) sf_error("No dz in vel");
    if(! sf_histfloat(vel,"o1",&oz)) sf_error("No oz in vel");
    if(! sf_histint(  vel,"n2",&nx)) sf_error("No nx in vel");  
    if(! sf_histfloat(vel,"d2",&dx)) sf_error("No dx in vel");
    if(! sf_histfloat(vel,"o2",&ox)) sf_error("No ox in vel");
    
    idx = 1/(dx*dx);
    idz = 1/(dz*dz);
        
    /* 
     * read velocity, density, reflectivity
     */
    vv=sf_floatalloc2(nz,nx);
    sf_floatread(vv[0],nz*nx,vel);
    
    dd=sf_floatalloc2(nz,nx);
    sf_floatread(dd[0],nz*nx,den);
    
    rr=sf_floatalloc2(nz,nx);
    sf_floatread(rr[0],nz*nx,ref);
    
    
    /* 
     * setup output header
     */
    sf_putint  (out,"n1",nz);
    sf_putfloat(out,"o1",dz);
    sf_putfloat(out,"d1",oz);
    sf_putint  (out,"n2",nx);
    sf_putfloat(out,"o2",dx);
    sf_putfloat(out,"d2",ox);
    sf_putint  (out,"n3",nt);
    sf_putfloat(out,"o3",dt);
    sf_putfloat(out,"d3",ot);

    

    /* 
     * allocate temporary arrays
     */

    up=sf_floatalloc2(nz,nx);
    uo=sf_floatalloc2(nz,nx);
    um=sf_floatalloc2(nz,nx);
    ud=sf_floatalloc2(nz,nx);
    
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ud[ix][iz]=0;
	}
    }
    
    fprintf(stderr,"%d %d \n",nz,nx);
    fprintf(stderr,"%f %f %f\n", c0,c1,c2);
    fprintf(stderr,"%f %f   \n", idx,idz);
    
    
    /* 
     *  MAIN LOOP
     */
    for (it=0; it<nt; it++) {
	fprintf(stderr,"%d\t",it);
	
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
	
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ud[ix][iz] *= vv[ix][iz]*vv[ix][iz];
	    }
	}
	
	/* exploding reflectors */
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		ud[ix][iz] -= wvl[it] * rr[ix][iz];
	    }
	}
	fprintf(stderr,"%g\t",wvl[it]);
	fprintf(stderr,"%g\n",ud[10][10]);
	
	/* time step */
	for (iz=0; iz<nz; iz++) {
	    for (ix=0; ix<nx; ix++) {
		up[ix][iz] = 2*uo[ix][iz] 
		    - um[ix][iz] 
		    + ud[ix][iz] * dt*dt;
		
		um[ix][iz] = uo[ix][iz];
		uo[ix][iz] = up[ix][iz];
	    }
	}
	
	/* write wavefield to output */
	sf_floatwrite(ud[0],nz*nx,out);
    }
    
    exit (0);
}
