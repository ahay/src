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
    off_t eseek;
    long int n2buf,i2buf,n1off;
    long int i1,i2,ii;

    float      **rin=NULL; /* data in (i1,n2buf) */
    sf_complex **cin=NULL;
    
    float       *rou=NULL; /* data out   (n2buf) */
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

    /* flag complex file */
    cmpl= (SF_COMPLEX == sf_gettype(Fin)?1:0);
    if(verb) sf_warning("cmplx=%d",cmpl);

    /* usable memory (Mb) */
    if (!sf_getfloat("memsize",&memsize)) memsize=1000.0;
    if(verb) sf_warning("memsize=%g",memsize);
     
    /*------------------------------------------------------------*/
    /* n2buf - how much of n2 can we load in memory? */
    n2buf=1;
    if(cmpl) while(n2buf/1024.*sf_n(a1)/1024.*SF_COMPLEX< memsize) n2buf++;
    else     while(n2buf/1024.*sf_n(a1)/1024.*SF_FLOAT  < memsize) n2buf++;
    
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
    /* allocate data arrays */

    if(verb) sf_warning("allocate arrays");
	
    if(cmpl) {
	cin=sf_complexalloc2(sf_n(a1),n2buf);
	cou=sf_complexalloc          (n2buf);
    } else {
	rin=sf_floatalloc2  (sf_n(a1),n2buf);
	rou=sf_floatalloc            (n2buf);
    }

    if(verb) sf_warning("OK");

    /* feedback index */
    ii = (long int)(sf_n(a1)/10.);
    
    /*------------------------------------------------------------*/
    n1off=0;
    for (i2=sf_n(a2); i2 > 0; i2 -= n2buf) {
        if (n2buf > i2) n2buf=i2;

	if(cmpl) {

	    /* read from input file */
	    sf_complexread  (cin[0],sf_n(a1)*n2buf,Fin);
	    
	    for(i1=0; i1<sf_n(a1); i1++) {
		
		/* copy data to output array - use pointers? */
		for(i2buf=0; i2buf<n2buf; i2buf++)
		    cou[i2buf] = cin[i2buf][i1];
		
		/* write to output file */
		eseek = (i1*sf_n(a2)+n1off)*sizeof(sf_complex);
		sf_seek(Fou,eseek,SEEK_SET);
		sf_complexwrite(cou,n2buf,Fou);
	    }

	} else {

	    /* read from input file */
	    sf_floatread  (rin[0],sf_n(a1)*n2buf,Fin);
	    
	    for(i1=0; i1<sf_n(a1); i1++) {
		
		/* copy data to output array - use pointers? */
		for(i2buf=0; i2buf<n2buf; i2buf++)
		    rou[i2buf] = rin[i2buf][i1];
		
		/* write to output file */
		eseek = (i1*sf_n(a2)+n1off)*sizeof(float);
		sf_seek(Fou,eseek,SEEK_SET);
		sf_floatwrite(rou,n2buf,Fou);

		if(verb && i1%ii==0) fprintf(stderr,".");
	    }
	}
	
	n1off+=n2buf;
	
	if(verb) fprintf(stderr," %5.1f%% complete\n",100.0*n1off/sf_n(a2));
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
