/* Convert CMPs to shots for regular 2-D geometry.

Takes: < cmps.rsf > shots.rsf

Assumes dy = dh!
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
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&oy)) sf_error("No o3= in input");

    ds = dy;
    os = oy + oh;
    ns = ny + nh - 1;

    sf_putint(out,"n3",ns);
    sf_putfloat(out,"d3",ds);
    sf_putfloat(out,"o3",os);

    trace = sf_floatalloc(nt);
    zero = sf_floatalloc(nt);
    for (it=0; it < nt; it++) {
	zero[it] = 0.;
    }

    tsize = nt*sizeof(float);

    sf_unpipe(in,ny*nh*tsize);
    pos = sf_tell(in);

    for (is=0; is < ns; is++) {
	sf_warning("shot %d of %d",is+1,ns);
	for (ih=0; ih < nh; ih++) {
	    iy = is - ih;
	    if (iy >= 0 && iy < ny) {
		sf_seek(in,pos+(iy*nh+ih)*tsize,SEEK_SET);
		sf_read(trace,sizeof(float),nt,in);
		sf_write(trace,sizeof(float),nt,out);
	    } else {
		sf_write(zero,sizeof(float),nt,out);
	    }
	}
    }

    exit(0);
}

