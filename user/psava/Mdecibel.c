/* Decibel */
#include <rsf.h>
#include <math.h>

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    bool verb;     /* verbosity flag */
    bool inv;

    sf_file  Fin, Fou; /* I/O files */

    sf_axis a1; /* cube axes */
    int     n,nbuf,ibuf;

    float aref;
    float *trace;

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getfloat("aref",&aref)) aref=1.0;   /* reference amplitude */
    if(! sf_getbool("inv",&inv)) inv=false; /* inverse transform */

    Fin = sf_input ("in" );
    Fou = sf_output("out");

    a1=sf_iaxa(Fin,1); if(verb) sf_raxa(a1); 

    nbuf = sf_bufsiz(Fin)/sizeof(float);
    trace = sf_floatalloc (nbuf);

    for (n = sf_filesize(Fin); n > 0; n -= nbuf) {
        if (nbuf > n) nbuf=n;

        sf_floatread(trace,nbuf,Fin);

	if(inv) {
	    for(ibuf=0; ibuf<nbuf; ibuf++)
		trace[ibuf] = aref * pow( 10,0.05*trace[ibuf]);
	} else {
	    for(ibuf=0; ibuf<nbuf; ibuf++)
		trace[ibuf] = 20. * log10( fabs(trace[ibuf]/aref) );
	}

        sf_floatwrite(trace,nbuf,Fou);

    }

    free(trace);

}
