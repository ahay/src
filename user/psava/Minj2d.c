/* inject/extract in/from wavefield */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fdutil.h"

int main(int argc, char* argv[])
{
    bool verb,adj; 

    /* I/O files */
    sf_file Ftrc=NULL; /* traces       */
    sf_file Fcoo=NULL; /* coordinates  */
    sf_file Fwfl=NULL; /* wavefield    */

    /* cube axes */
    sf_axis at,az,ax,aa;
    sf_axis aco;

    /* I/O arrays */
    float  *wco=NULL;  /* traces   */
    pt2d   *coo=NULL;  /* coordinates   */
    lint2d  cow;       /* weights/indices */
    float **wfl=NULL;  /* wavefield   */

    fdm2d    fdm=NULL;
    int nco;
    int   iz,ix,it;
    int   nz,nx,nt;
    float dz,dx;
    float oz,ox;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */
    if(! sf_getbool("adj", &adj))   adj=false; /* adjoint flag */

    /*------------------------------------------------------------*/
    /* I/O files */
    
    Fcoo = sf_input ("coo"); /* coordinates */
    sf_oaxa(Fcoo,aco,1);
    nco = sf_n(aco);
    coo = (pt2d*) sf_alloc(nco,sizeof(*coo)); 
    pt2dread1(Fcoo,coo,nco,2); /* read (x,z) coordinates */

    if(adj) {
	Fwfl = sf_input("wfl"); /* wavefield */
	Ftrc = sf_output("in" ); /* traces   */

	az = sf_iaxa(Fwfl,1); sf_setlabel(az,"z");
	ax = sf_iaxa(Fwfl,2); sf_setlabel(ax,"x");
	at = sf_iaxa(Fwfl,3); sf_setlabel(at,"t");

	aa = sf_maxa(1,0,1);
	sf_oaxa(Ftrc,aco,1);
	sf_oaxa(Ftrc,at,2);
	sf_oaxa(Ftrc,aa,3);

    } else {
	Ftrc = sf_input ("in" ); /* traces   */
	Fwfl = sf_output("wfl"); /* wavefield */

	at = sf_iaxa(Ftrc,2); sf_setlabel(at,"t");

	if(!sf_getint  ("nz",&nz)) nz=1;
	if(!sf_getfloat("oz",&oz)) oz=0.0;
	if(!sf_getfloat("dz",&dz)) dz=1.0;
	sf_setlabel(az,"z");
	az = sf_maxa(nz,oz,dz);

	if(!sf_getint  ("nx",&nx)) nx=1; 
	if(!sf_getfloat("ox",&ox)) ox=0.0;
	if(!sf_getfloat("dx",&dx)) dx=1.0;
	sf_setlabel(ax,"x");
	ax = sf_maxa(nx,ox,dx);

	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Ftrc,at,3);
    }
    
    if(verb) {
	sf_raxa(az);
	sf_raxa(ax);
	sf_raxa(at);
    }


    /* allocate wavefield arrays */
    wco = sf_floatalloc(nco);
    wfl = sf_floatalloc2(sf_n(az),sf_n(ax));
    for    (ix=0; ix<sf_n(ax); ix++) {
	for(iz=0; iz<sf_n(az); iz++) {
	    wfl[ix][iz]=0;
	}
    }

    /* interpolation coefficients */
    fdm = fdutil_init(verb,'n',az,ax,0,1);
    cow = lint2d_make(nco,coo,fdm);
 
    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
	
	if(adj) {
	    sf_floatread(wfl[0],sf_n(az)*sf_n(ax),Fwfl);
	    lint2d_extract(wfl,wco,cow);
	    sf_floatwrite(wco,nco,Ftrc);
	} else {
	    sf_floatread(wco,nco,Ftrc);
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
    /* CLOSE FILES AND EXIT */
    if (Ftrc!=NULL) sf_fileclose(Ftrc); 
    if (Fwfl!=NULL) sf_fileclose(Fwfl);
    if (Fcoo!=NULL) sf_fileclose(Fcoo);

    exit (0);
}
