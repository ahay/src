#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2, it;
    float t, *data, *slope, *next;
    sf_file trace, dip, out;

    sf_init(argc,argv);
    trace = sf_input("trace");
    dip = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in dip");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in dip");
   
    data = sf_floatalloc(n1);
    next = sf_floatalloc(n1);
    slope = sf_floatalloc(n1);

    sf_floatread(next,n1,trace);
    sf_floatwrite(next,n1,out);
    
    for (i2=1; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    data[i1] = next[i1];
	}
	sf_floatread(slope,n1,dip);
	for (i1=0; i1 < n1; i1++) {
	    t = i1-slope[i1];
	    it = t; t -= it;
	    if (it <= 0) {
		next[i1] = data[0];
	    } else if (it >= n1-1) {
		next[i1] = data[n1-1];
	    } else {
		next[i1] = t*data[it+1] + (1.-t)*data[it];
	    }
	}
	sf_floatwrite(next,n1,out);
    }


    exit(0);
}

/* 	$Id$	 */
