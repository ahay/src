/* Convert CMPs to shots for regular 2-D geometry.

Takes: < cmps.rsf > shots.rsf

Assumes ds=dh.
*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt,ns, ny,nh, iy,ih,is,it;
    long pos, tsize;
    float ds, dy,dh, os, oy,oh;
    float *trace, *zero;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ns)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&ds)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&os)) sf_error("No o3= in input");

    dy = ds;
    oy = os + oh + (nh-1)*dh;
    ny = ns - nh;

    sf_putint(out,"n3",ny);
    sf_putfloat(out,"d3",dy);
    sf_putfloat(out,"o3",oy);

    trace = sf_floatalloc(nt);
    zero = sf_floatalloc(nt);
    for (it=0; it < nt; it++) {
	zero[it] = 0.;
    }

    tsize = nt*sizeof(float);

    sf_unpipe(in,ns*nh*tsize);
    pos = sf_tell(in);

    for (iy=0; iy < ny; iy++) {
	for (ih=0; ih < nh; ih++) {
	    is = iy - ih + nh - 1;
	    if (is >= 0 && is < ns) {
		sf_seek(in,pos+(is*nh+ih)*tsize,SEEK_SET);
		sf_read(trace,sizeof(float),nt,in);
		sf_write(trace,sizeof(float),nt,out);
	    } else {
		sf_write(zero,sizeof(float),nt,out);
	    }
	}
    }

    exit(0);
}

