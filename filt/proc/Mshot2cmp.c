/* Convert shots to CMPs for regular 2-D geometry.

Takes: < cmps.rsf > shots.rsf
*/
#include <string.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt,ns, ny,nh, iy,ih,is, type, esize;
    long pos;
    float ds, dy,dh, os, oy,oh;
    char *trace, *zero;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ns)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&ds)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&os)) sf_error("No o3= in input");

    type = 0.5 + ds/dh;

    dy = dh;
    oy = os + oh;
    ny = ns*type + nh - 1;

    sf_putint(out,"n2",(nh+type-1)/type);
    sf_putfloat(out,"d2",type*dh);

    sf_putint(out,"n3",ny);
    sf_putfloat(out,"d3",dy);
    sf_putfloat(out,"o3",oy);

    if (!sf_histint(in,"esize",&esize)) {
	esize=4;
    } else if (0>=esize) {
	sf_error("wrong esize=%d",esize);
    }
    nt *= esize;

    trace = sf_charalloc(nt);
    zero = sf_charalloc(nt);
    memset(zero,0,nt);

    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");
    
    sf_unpipe(in,(long) ns*nh*nt);
    pos = sf_tell(in);

    for (iy=0; iy < ny; iy++) {
	for (ih=iy%type; ih < nh; ih += type) {
	    is = (iy - ih)/type;
	    if (is >= 0 && is < ns) {
		sf_seek(in,pos+(is*nh+ih)*nt,SEEK_SET);
		sf_read(trace,sizeof(char),nt,in);
		sf_write(trace,sizeof(char),nt,out);
	    } else {
		sf_write(zero,sizeof(char),nt,out);
	    }
	}
    }

    exit(0);
}

