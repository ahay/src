/* Window a data set based on a header mask.

Takes: < input.rsf mask=mask.rsf > windowed.rsf

The input data is 2-D collection of traces n1xn2,
mask is a 1-D integer array n2, windowed is n1xm2,
where m2 is the number of nonzero elements in mask.
*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, j2, i2, esize, *mask;
    long pos;
    char *trace;
    sf_file in, head, out;

    sf_init (argc,argv);
 
    head = sf_input("mask");
    if (SF_INT != sf_gettype(head))
	sf_error("Need integer mask");
    n2 = sf_filesize(head);
 
    mask = sf_intalloc(n2);
    
    sf_read(mask,sizeof(int),n2,head);
    sf_fileclose(head);

    in = sf_input ("in");
    out = sf_output ("out");
 
    if (!sf_histint(in,"n1",&n1)) n1=1;
    if (!sf_histint(in,"esize",&esize)) esize=4;
    n1 *= esize;

    trace = sf_charalloc(n1);

    if (!sf_histint(in,"n2",&j2) || j2 != n2 || sf_leftsize(in,1) != n2)
	sf_error("Wrong input dimensions, need %d by %d",n1,n2);
    
    for (j2=i2=0; i2 < n2; i2++) {
	if (mask[i2]) j2++;
    }
    sf_putint(out,"n2",j2);

    sf_unpipe(in,n1*n2);
    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");

    pos = sf_tell(in);
    for (i2=0; i2<n2; i2++) {
	if (mask[i2]) {
	    sf_seek(in,pos+i2*n1,SEEK_SET);
	    sf_read(trace,1,n1,in);
	    sf_write(trace,1,n1,out);
	}
    }

    exit(0);
}
    
