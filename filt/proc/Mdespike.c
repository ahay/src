/* Remove spikes in by sliding 1-D medians. */

#include <rsf.h>

#include "quantile.h"

int main(int argc, char* argv[]) 
{
    int n1, wide, i1, shift, i, i2, n2;
    float *data, *signal, *win;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc(n1);
    signal = sf_floatalloc(n1);

    if (!sf_getint("wide",&wide)) wide=7;
    /* sliding window width */

    win = sf_floatalloc(wide);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);
	
	for (i1=0; i1 < n1; i1++) {
	    shift = i1-wide/2;
	    if (shift+wide-1 >= n1) {
		shift = n1-wide;
	    } else if (shift < 0) {
		shift = 0;
	    }
	    for (i=0; i < wide; i++) {
		win[i] = data[shift+i];
	    }
	    signal[i1] = quantile(1+wide/2,wide,win);
	}
	
	sf_floatwrite(signal,n1,out);
    }	

    sf_close();
    exit(0);
}
