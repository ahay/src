/* inject/extract in/from 3D wavefield */
#include <rsf.h>
#include "fdutil.h"

int main(int argc, char* argv[])
{
    bool verb,adj; 

    /* I/O files */
    sf_file Ftrc=NULL; /* traces       */
    sf_file Fcoo=NULL; /* coordinates  */
    sf_file Fwfl=NULL; /* wavefield    */

    /* cube axes */
    sf_axis at,az,ax,ay,aa,ac;

    /* I/O arrays */
    float   *wco=NULL;  /* traces   */
    pt3d    *coo=NULL;  /* coordinates   */
    lint3d   cow;       /* weights/indices */
    float ***wfl=NULL;  /* wavefield   */

    fdm3d fdm=NULL;
    int   iz,iy,ix,it;
    int   nz,ny,nx;
    float dz,dy,dx;
    float oz,oy,ox;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("adj", &adj))   adj=false; /* adjoint flag */

    /*------------------------------------------------------------*/
    /* setup I/O */
    Fcoo = sf_input("coo"); /* coordinates */
    ac = sf_iaxa(Fcoo,3); sf_setlabel(ac,"c"); sf_setunit(ac,"");
    coo = (pt3d*) sf_alloc(sf_n(ac),sizeof(*coo)); 
    pt3dread1(Fcoo,coo,sf_n(ac),3); /* read (x,y,z) coordinates */

    if(adj) {
	Fwfl = sf_input ("in");  /* wavefield */
	Ftrc = sf_output("out"); /* traces   */

	az = sf_iaxa(Fwfl,1); sf_setlabel(az,"z");
	ax = sf_iaxa(Fwfl,2); sf_setlabel(ax,"x");
	ay = sf_iaxa(Fwfl,3); sf_setlabel(ay,"y");
	at = sf_iaxa(Fwfl,4); sf_setlabel(at,"t");

	aa = sf_maxa(1,0,1);
	sf_oaxa(Ftrc,ac,1);
	sf_oaxa(Ftrc,at,2);
	sf_oaxa(Ftrc,aa,3);
	sf_oaxa(Ftrc,aa,4);
    } else {
	Ftrc = sf_input ("in" ); /* traces   */
	Fwfl = sf_output("out"); /* wavefield */

	at = sf_iaxa(Ftrc,2); sf_setlabel(at,"t");

	if(!sf_getint  ("nz",&nz)) nz=1;
	if(!sf_getfloat("oz",&oz)) oz=0.0;
	if(!sf_getfloat("dz",&dz)) dz=1.0;
	az = sf_maxa(nz,oz,dz);
	sf_setlabel(az,"z");

	if(!sf_getint  ("ny",&ny)) ny=1; 
	if(!sf_getfloat("oy",&oy)) oy=0.0;
	if(!sf_getfloat("dy",&dy)) dy=1.0;
	ay = sf_maxa(ny,oy,dy);
	sf_setlabel(ay,"y");

	if(!sf_getint  ("nx",&nx)) nx=1; 
	if(!sf_getfloat("ox",&ox)) ox=0.0;
	if(!sf_getfloat("dx",&dx)) dx=1.0;
	ax = sf_maxa(nx,ox,dx);
	sf_setlabel(ax,"x");

	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Fwfl,ay,3);
	sf_oaxa(Fwfl,at,4);
    }
    
    if(verb) {
	sf_raxa(az);
	sf_raxa(ay);
	sf_raxa(ax);
	sf_raxa(at);
	sf_raxa(ac);	
    }

    /* allocate wavefield arrays */
    wco = sf_floatalloc (sf_n(ac));
    wfl = sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay));

    /* interpolation coefficients */
    fdm = fdutil3d_init(verb,'n',az,ax,ay,0,1);
    cow = lint3d_make(sf_n(ac),coo,fdm);

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<sf_n(at); it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	if(adj) {
	    sf_floatread(wfl[0][0],sf_n(az)*sf_n(ay)*sf_n(ax),Fwfl);

	    lint3d_extract(wfl,wco,cow);

	    sf_floatwrite(wco,sf_n(ac),Ftrc);
	} else {
	    sf_floatread(wco,sf_n(ac),Ftrc);

	    for    (iy=0; iy<sf_n(ay); iy++)
		for    (ix=0; ix<sf_n(ax); ix++)
		    for(iz=0; iz<sf_n(az); iz++)
			wfl[iy][ix][iz]=0;
	    lint3d_inject(wfl,wco,cow);

	    sf_floatwrite(wfl[0][0],sf_n(az)*sf_n(ay)*sf_n(ax),Fwfl);
	}

    }	/* end time loop */
    if(verb) fprintf(stderr,"\n");
	
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(**wfl); free(*wfl); free(wfl);
    free(wco);
    free(coo);
    
    /*------------------------------------------------------------*/ 
    /* close files */
    if (Ftrc!=NULL) sf_fileclose(Ftrc); 
    if (Fwfl!=NULL) sf_fileclose(Fwfl);
    if (Fcoo!=NULL) sf_fileclose(Fcoo);

    exit (0);
}
