#include <rsf.h>
#include "wefic.h"

int main (int argc, char *argv[])
{
    bool  verb;           /* verbosity */
    char *itype;          /* imaging type (zero-time, finite-time) */
    axa az,ax,ay,aw,ah;

    sf_file Fus; /* source   wavefield D(x,y,z,w) */
    sf_file Fur; /* receiver wavefield U(x,y,z,w) */
    sf_file Fi;  /*    image           R(x,y,z,h) */

    fslice imag,sdat,rdat;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (!sf_getbool("verb",&verb)) verb =  true; /* verbosity flag */    
    if (NULL == (itype = sf_getstring("itype"))) itype = "z";

    Fus = sf_input ( "in");
    Fur = sf_input ("rwf");
    Fi  = sf_output("out");

    iaxa(Fus,&ax,1); ax.l="x"; oaxa(Fi,&ax,1);
    iaxa(Fus,&ay,2); ay.l="y"; oaxa(Fi,&ay,2);    
    iaxa(Fus,&az,3); az.l="z"; oaxa(Fi,&az,3);
    iaxa(Fus,&aw,4); aw.l="w";
	    
    sdat = fslice_init( ax.n*ay.n*az.n*aw.n,1,sizeof(float complex));
    rdat = fslice_init( ax.n*ay.n*az.n*aw.n,1,sizeof(float complex));

    switch(itype[0]) {
	case 'p':
	    sf_warning("Full offset IC");
	    sf_settype(Fi,SF_FLOAT);	    

	    if(!sf_getint  ("nh",&ah.n)) ah.n=1;
	    if(!sf_getfloat("oh",&ah.o)) ah.o=0;
	    if(!sf_getfloat("dh",&ah.d)) ah.d=0.1;
	    ah.l="ht";

	    oaxa(Fi,&ah,4);

	    imag = fslice_init( ax.n*ay.n*az.n,ah.n,sizeof(float));
	    break;
	case 'z':
	default:
	    sf_warning("Zero offset IC");
	    sf_settype(Fi,SF_FLOAT);

	    ah.n=1; ah.o=0; ah.d=1; ah.l=" ";
	    oaxa(Fi,&ah,4);

	    imag = fslice_init( ax.n*ay.n*az.n,1,sizeof(float));
	    break;
    }

    fslice_load(Fus,sdat,SF_COMPLEX);
    fslice_load(Fur,rdat,SF_COMPLEX);
    /*------------------------------------------------------------*/
    wefic_init(verb,ax,ay,az,aw);

    switch(itype[0]) {
	case 'p':
	    hhfic_init(ah);
	    hhfic(sdat,rdat,imag);
	    hhfic_close();
	    break;
	case 'z':
	default:
	    zofic(sdat,rdat,imag);
	    break;
    }
    wefic_close();
    /*------------------------------------------------------------*/    

    /* slice management (temp files) */
    fslice_dump(Fi,imag,SF_FLOAT);

    fslice_close(sdat);
    fslice_close(rdat);
    fslice_close(imag);

    exit (0);
}
