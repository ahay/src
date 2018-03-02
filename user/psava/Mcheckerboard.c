/* make a 2D/3D checkerboard model */
#include <rsf.h>

int main(int argc, char* argv[])
{
    int   nx, ny, nz;
    float ox, oy, oz;
    float dx, dy, dz;
    int   ix, iy, iz;
    bool  cx, cy, cz;
    int N;
    
    float ***ccc;

    sf_file Fo;

    sf_init (argc, argv);
    

    /*------------------------------------------------------------*/
    /* get dimensions */
    if(! sf_getint("nx",&nx)) nx=1;
    if(! sf_getint("ny",&ny)) ny=1;
    if(! sf_getint("nz",&nz)) nz=1;
    
    if(! sf_getfloat("ox",&ox)) ox=0.0;
    if(! sf_getfloat("oy",&oy)) oy=0.0;
    if(! sf_getfloat("oz",&oz)) oz=0.0;

    if(! sf_getfloat("dx",&dx)) dx=1.0;
    if(! sf_getfloat("dy",&dy)) dy=1.0;
    if(! sf_getfloat("dz",&dz)) dz=1.0;

    /* get checkerboard size */
    if(! sf_getint("N",&N)) N=1;
    
    /*------------------------------------------------------------*/
    Fo = sf_output("out");

    /* make output header */
    sf_putint  (Fo,"n1",nz);
    sf_putfloat(Fo,"d1",dz);
    sf_putfloat(Fo,"o1",oz);

    sf_putint  (Fo,"n2",nx);
    sf_putfloat(Fo,"d2",dx);
    sf_putfloat(Fo,"o2",ox);

    sf_putint  (Fo,"n3",ny);
    sf_putfloat(Fo,"d3",dy);
    sf_putfloat(Fo,"o3",oy);
    /*------------------------------------------------------------*/

    ccc = sf_floatalloc3( nz, nx, ny);
    
    /*------------------------------------------------------------*/

    cy = false;
    for(iy=0; iy<ny; iy++) {
	if( iy%N ==0 ) cy = !cy;

	cx = true;
	for(ix=0; ix<nx; ix++) {
	    if( ix%N ==0 ) cx = !cx;

	    cz = true;
	    for(iz=0; iz<nz; iz++) {
		if( iz%N ==0 ) cz = !cz;

		if(  ( cy && (  cx &&  cz) ) ||
		     ( cy && ( !cx && !cz) ) ||
		     (!cy && (  cx && !cz) ) ||
		     (!cy && ( !cx &&  cz) ) )
		    ccc[iy][ix][iz]=1;
			
	    }
	}
    }
    
    /*------------------------------------------------------------*/
    sf_floatwrite(ccc[0][0],nz*nx*ny,Fo);
    
    exit (0);
}

