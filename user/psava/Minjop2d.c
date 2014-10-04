/* inject/extract in/from 2D wavefield */
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
    sf_axis at,az,ax,aa,ac;

    /* I/O arrays */
    float  *wco=NULL;  /* traces   */
    pt2d   *coo=NULL;  /* coordinates   */
    lint2d  cow;       /* weights/indices */
    float **wfl=NULL;  /* wavefield   */

    fdm2d fdm=NULL;
    int   iz,ix,it;
    int   nz,nx;
    float dz,dx;
    float oz,ox;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("adj", &adj))   adj=false; /* adjoint flag */

    /*------------------------------------------------------------*/
    /* setup I/O */
    Fcoo = sf_input ("coo"); /* coordinates */
    ac = sf_iaxa(Fcoo,2);
    sf_setlabel(ac,"c");
    sf_setunit(ac,"");
    coo = (pt2d*) sf_alloc(sf_n(ac),sizeof(*coo)); 
    pt2dread1(Fcoo,coo,sf_n(ac),2); /* read (x,z) coordinates */

    if(adj) {
	Fwfl = sf_input ("in");  /* wavefield */
	Ftrc = sf_output("out"); /* traces   */

	az = sf_iaxa(Fwfl,1); sf_setlabel(az,"z");
	ax = sf_iaxa(Fwfl,2); sf_setlabel(ax,"x");
	at = sf_iaxa(Fwfl,3); sf_setlabel(at,"t");

	aa = sf_maxa(1,0,1);
	sf_oaxa(Ftrc,ac,1);
	sf_oaxa(Ftrc,at,2);
	sf_oaxa(Ftrc,aa,3);
    } else {
	Ftrc = sf_input ("in" ); /* traces   */
	Fwfl = sf_output("out"); /* wavefield */

	at = sf_iaxa(Ftrc,2); sf_setlabel(at,"t");

	if(!sf_getint  ("nz",&nz)) nz=1;
	if(!sf_getfloat("oz",&oz)) oz=0.0;
	if(!sf_getfloat("dz",&dz)) dz=1.0;
	az = sf_maxa(nz,oz,dz);
	sf_setlabel(az,"z");

	if(!sf_getint  ("nx",&nx)) nx=1; 
	if(!sf_getfloat("ox",&ox)) ox=0.0;
	if(!sf_getfloat("dx",&dx)) dx=1.0;
	ax = sf_maxa(nx,ox,dx);
	sf_setlabel(ax,"x");

	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Fwfl,at,3);
    }
    
    if(verb) {
	sf_raxa(az);
	sf_raxa(ax);
	sf_raxa(at);
	sf_raxa(ac);	
    }

    /* allocate wavefield arrays */
    wco = sf_floatalloc (sf_n(ac));
    wfl = sf_floatalloc2(sf_n(az),sf_n(ax));

    /* interpolation coefficients */
    fdm = fdutil_init(verb,'n',az,ax,0,1);
    cow = lint2d_make(sf_n(ac),coo,fdm);

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<sf_n(at); it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	if(adj) {
	    sf_floatread(wfl[0],sf_n(az)*sf_n(ax),Fwfl);

	    lint2d_extract(wfl,wco,cow);

	    sf_floatwrite(wco,sf_n(ac),Ftrc);
	} else {
	    sf_floatread(wco,sf_n(ac),Ftrc);

	    for    (ix=0; ix<sf_n(ax); ix++)
		for(iz=0; iz<sf_n(az); iz++)
		    wfl[ix][iz]=0;
	    lint2d_inject(wfl,wco,cow);

	    sf_floatwrite(wfl[0],sf_n(az)*sf_n(ax),Fwfl);
	}

    }	/* end time loop */
    if(verb) fprintf(stderr,"\n");
	
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*wfl); free(wfl);
    free(wco);
    free(coo);
    
    /*------------------------------------------------------------*/ 
    /* close files */
    if (Ftrc!=NULL) sf_fileclose(Ftrc); 
    if (Fwfl!=NULL) sf_fileclose(Fwfl);
    if (Fcoo!=NULL) sf_fileclose(Fcoo);

    exit (0);
}
