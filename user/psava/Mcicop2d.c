/* Conventional IC 2D */
/* 
 * ----------------------------------------
 * FORWARD OP: opr * wfl = img
 *     inputs:   z- x -t   opr 
 *               z- x- t   wfl
 *     output:   z- x      img
 * ----------------------------------------
 * ADJOINT OP: opr * img = wfl
 *     inputs:   z- x- t   opr
 *               z- x      img
 *     output:   z- x- t   wfl
 * ----------------------------------------
 */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define CICLOOP(a)				\
    for    (ix=0;ix<nx;ix++){                   \
        for(iz=0;iz<nz;iz++){                   \
            {a}                                 \
        }                                       \
    }

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    bool verb;     /* verbosity flag */
    bool adj;      /* adjoint operator flag */
    bool wflcausal, oprcausal; /* causal wfl?, opr? */

    sf_file  Fopr,      Fwfl,       Fimg; /* I/O files */
    float   **opr=NULL,**wfl=NULL,**img=NULL;

    sf_axis az,ax,at,aa; /* cube axes */
    int     iz,ix,it;
    int     nz,nx,nt;

    float scale; /* time summation scaling */
    int nslice;  /* wavefield slice size */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool(    "verb",&verb    ))        verb=false; /* verbosity flag */
    if(! sf_getbool(     "adj",&adj     ))         adj=false; /* adjoint flag */
    if(! sf_getbool("wflcausal",&wflcausal)) wflcausal=false; /* causal wfl? */
    if(! sf_getbool("oprcausal",&oprcausal)) oprcausal=false; /* causal opr? */

    /*------------------------------------------------------------*/
    Fopr = sf_input ("opr" ); /*   operator */
    az=sf_iaxa(Fopr,1); if(verb) sf_raxa(az); nz = sf_n(az);
    ax=sf_iaxa(Fopr,2); if(verb) sf_raxa(ax); nx = sf_n(ax);
    at=sf_iaxa(Fopr,3); if(verb) sf_raxa(at); nt = sf_n(at);
    aa=sf_maxa(1,0,1); sf_setlabel(aa,"");  sf_setunit (aa,""); 
    
    scale = 1./nt;               /* time summation scaling */
    nslice = nz*nx*sizeof(float); /* wavefield slice */

    /*------------------------------------------------------------*/
    /* setup output */
    if(adj) {
	Fimg = sf_input ( "in"); /*  read img */

	Fwfl = sf_output("out"); /* write wfl */
	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Fwfl,at,3);
    } else {
	Fwfl = sf_input ( "in"); /*  read wfl */

	Fimg = sf_output("out"); /* write img */
	sf_oaxa(Fimg,az,1);
	sf_oaxa(Fimg,ax,2);
	sf_oaxa(Fimg,aa,3);
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate arrays */
    img = sf_floatalloc2(nz,nx); 
    opr = sf_floatalloc2(nz,nx);
    wfl = sf_floatalloc2(nz,nx);
    /*------------------------------------------------------------*/

    if(adj) { /* ADJOINT OPERATOR */
	for(it=0;it<nt;it++) sf_floatwrite(wfl[0],nz*nx,Fwfl); /* reserve wfl */ 
	sf_seek(Fwfl,0,SEEK_SET);                              /* seek back */

	sf_floatread(img[0],nz*nx,Fimg);   /*  read img */
	CICLOOP( img[ix][iz] *=scale; );   /* scale img */
	
	for(it=0;it<nt;it++){ if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	   	   
	    if( !oprcausal ) sf_seek(Fopr,(off_t)(nt-1-it)*nslice,SEEK_SET);
	    sf_floatread(opr[0],nz*nx,Fopr); /* read opr */

#ifdef _OPENMP	    
#pragma omp parallel for schedule(dynamic)				\
    private(ix,iz)							\
    shared (wfl,opr,img,nx,nz)
#endif
	    CICLOOP( wfl[ix][iz] = opr[ix][iz]*img[ix][iz]; ); /* CIC ADJ */

	    if( !wflcausal ) sf_seek(Fwfl,(off_t)(nt-1-it)*nslice,SEEK_SET);
	    sf_floatwrite(wfl[0],nz*nx,Fwfl); /* write wfl */
	    
	} /* it */
	if(verb) fprintf(stderr,"\n");

    } else { /* FORWARD OPERATOR */
	CICLOOP( img[ix][iz]=0.; ); /* zero img */

	sf_floatwrite(img[0],nz*nx,Fimg); /* reserve img */
	sf_seek(Fimg,0,SEEK_SET);         /* seek back */

	for(it=0; it<nt; it++) { if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	   	    
	    if( !oprcausal ) sf_seek(Fopr,(off_t)(nt-1-it)*nslice,SEEK_SET);
	    sf_floatread(opr[0],nz*nx,Fopr); /* read opr */

	    if( !wflcausal ) sf_seek(Fwfl,(off_t)(nt-1-it)*nslice,SEEK_SET);
	    sf_floatread(wfl[0],nz*nx,Fwfl); /* read wfl */
	    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ix,iz)							\
    shared (wfl,opr,img,nx,nz)
#endif
	    CICLOOP( img[ix][iz] += opr[ix][iz]*wfl[ix][iz]; ); /* CIC FOR */

	} /* it */
	if(verb) fprintf(stderr,"\n");
	
	CICLOOP( img[ix][iz] *=scale; );  /* scale img */
	sf_floatwrite(img[0],nz*nx,Fimg); /* write img */

    } /* end if adj */

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*img); free(img);
    free(*opr); free(opr);
    free(*wfl); free(wfl);
    /*------------------------------------------------------------*/
    
    exit (0);
}
