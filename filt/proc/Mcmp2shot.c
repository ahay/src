/* Convert CMPs to shots for regular 2-D geometry.

Takes: < cmps.rsf > shots.rsf

*/
#include <string.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt,ns, ny,nh, iy,ih,is,it, type, esize;
    long pos;
    bool sign;
    float ds, dy,dh, os, oy,oh;
    char *trace, *zero;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&oy)) sf_error("No o3= in input");

    if (!sf_getbool("positive",&sign)) sign=true;
    /* initial offset orientation */

    type = 0.5 + dh/dy;

    ds = dh;
    os = sign? oy - oh - (nh-1)*dh: oy + oh;
    ns = (ny-1)/type + nh;

    sf_putint(out,"n2",type*nh);
    sf_putfloat(out,"d2",dh/type);

    sf_putint(out,"n3",ns);
    sf_putfloat(out,"d3",ds);
    sf_putfloat(out,"o3",os);

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

    sf_unpipe(in,(long) ny*nh*nt);
    pos = sf_tell(in);

    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < type; it++) {
		iy = sign? it + type*(is + ih - nh + 1): type*(is - ih) - it;
		if (iy >= 0 && iy < ny) {
		    sf_seek(in,pos+(iy*nh+ih)*nt,SEEK_SET);
		    sf_read(trace,sizeof(char),nt,in);
		    sf_write(trace,sizeof(char),nt,out);
		} else {
		    sf_write(zero,sizeof(char),nt,out);
		}
	    }
	}
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mcmp2shot.c,v 1.4 2004/03/22 05:43:24 fomels Exp $	 */
