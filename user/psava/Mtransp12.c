/* Transpose 1-2 */
#include <rsf.h>
#include <math.h>

int main(int argc, char* argv[])
{
    bool verb; /* verbosity flag */
    bool cmpl; /* complex I/O */

    sf_file Fin,Fou; 
    sf_axis a1,a2;

    float  memsize; /* in Mb */
    size_t memelem,nbuf,noff,eseek;

    long int ibuf,i1,i2;

    float      **rin=NULL; /* (i1,nbuf) */
    float       *rou=NULL; /*    (nbuf) */
    sf_complex **cin=NULL;
    sf_complex  *cou=NULL;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;  /* verbosity flag */

    /*------------------------------------------------------------*/
    Fin = sf_input ("in");  /*  input file */
    a1 = sf_iaxa(Fin,1); 
    a2 = sf_iaxa(Fin,2);
    
    Fou = sf_output("out"); /* output file */
    sf_oaxa(Fou,a2,1);   
    sf_oaxa(Fou,a1,2);

    /* set complex flag */
    cmpl= (SF_COMPLEX == sf_gettype(Fin)?1:0);
    sf_warning("cmplx=%d",cmpl);

    /*------------------------------------------------------------*/
    if(verb) sf_warning("reserve output");

    if(cmpl) {
	cou=sf_complexalloc (sf_n(a2));
	for(i2=0;i2<sf_n(a2);i2++) cou[i2]=0.;
	for(i1=0;i1<sf_n(a1);i1++) sf_complexwrite(cou,sf_n(a2),Fou);
	free(cou);
    } else {
	rou=sf_floatalloc (sf_n(a2));
	for(i2=0;i2<sf_n(a2);i2++) rou[i2]=0.;
	for(i1=0;i1<sf_n(a1);i1++) sf_floatwrite(rou,sf_n(a2),Fou);
	free(rou);
    }
    sf_seek(Fou,0,SEEK_SET); 

    if(verb) sf_warning("OK");

    /*------------------------------------------------------------*/
    /*    sf_warning("sizeof(size_t) = %d",sizeof(memelem));*/
    /*    sf_warning("sizeof(int) = %d",sizeof(ibuf));*/

    /* usable memory Mb */
    if (!sf_getfloat("memsize",&memsize)) memsize=1000.0;

    /* usable memory elements */
    if(cmpl) memelem = memsize/(SF_COMPLEX)*1024*1024;
    else     memelem = memsize/(SF_FLOAT)  *1024*1024;
    /*    sf_warning("memelem=%zd",memelem);*/

    /* nbuf - how much of n2 can load we load in memory? */
    nbuf=SF_MIN( floor( (1.0*memelem/sf_n(a1) )) , sf_n(a2) );
    if(verb) sf_warning("nbuf=%d",nbuf);
    
    /*------------------------------------------------------------*/
    /* allocate data arrays */
    if(cmpl) {
	cin=sf_complexalloc2(sf_n(a1),nbuf); /* data in */
	cou=sf_complexalloc          (nbuf); /* data out */
    } else {
	rin=sf_floatalloc2(sf_n(a1),nbuf); /* data in */
	rou=sf_floatalloc          (nbuf); /* data out */
    }

    /*------------------------------------------------------------*/
    noff=0;
    for (i2=sf_n(a2); i2 > 0; i2 -= nbuf) {
        if (nbuf > i2) nbuf=i2;

	if(cmpl) {

	    /* read from input file */
	    sf_complexread  (cin[0],sf_n(a1)*nbuf,Fin);
	    
	    for(i1=0; i1<sf_n(a1); i1++) {
		
		/* copy data to output array - use pointers? */
		for(ibuf=0; ibuf<nbuf; ibuf++)
		    cou[ibuf] = cin[ibuf][i1];
		
		/* write to output file */
		eseek = i1*sf_n(a2)+noff; /* seek elements */
		sf_seek(Fou,(off_t)(eseek*sizeof(sf_complex)),SEEK_SET);
		sf_complexwrite(cou,nbuf,Fou);
	    }

	} else {

	    /* read from input file */
	    sf_floatread  (rin[0],sf_n(a1)*nbuf,Fin);
	    
	    for(i1=0; i1<sf_n(a1); i1++) {
		
		/* copy data to output array - use pointers? */
		for(ibuf=0; ibuf<nbuf; ibuf++)
		    rou[ibuf] = rin[ibuf][i1];
		
		/* write to output file */
		eseek = i1*sf_n(a2)+noff; /* seek elements */
		sf_seek(Fou,(off_t)(eseek*sizeof(float)),SEEK_SET);
		sf_floatwrite(rou,nbuf,Fou);
	    }

	}

	noff+=nbuf;
	/*	sf_warning("noff=%d",noff);*/
	if(verb) sf_warning("%5.1f%% complete",100.0*noff/sf_n(a2));
    }

    /*------------------------------------------------------------*/
    if(cmpl) {
	;           free(cou);
	free(*cin); free(cin);
    } else {
	;           free(rou);
	free(*rin); free(rin);
    }

    exit (0);
}
