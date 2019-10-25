/* Apply running mean or median filter */

#include <rsf.h>

static float slow_median(int n, float* list)
/* find median by slow sorting, changes list */
{
    int k, k2;
    float item1, item2;
    
    for (k=0; k < n; k++) {
	item1 = list[k];

	/* assume everything up to k is sorted */
        for (k2=k; k2 > 0; k2--) {
            item2 = list[k2-1];
            if (item1 >= item2) break;
	    list[k2]   = item2;
	}
	list[k2] = item1;
    }
    
    return list[n/2];
}

int main(int argc, char* argv[]) 
{
    int w1, w2, nw, s1,s2, j1,j2, i1,i2,i3, n1,n2,n3;
    char *how;
    float **data, **signal, **win;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* get data dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2=");
    n3 = sf_leftsize(in,2);

    /* input and output */
    data = sf_floatalloc2(n1,n2);
    signal = sf_floatalloc2(n1,n2);

    if (!sf_getint("w1",&w1)) w1=5;
    if (!sf_getint("w2",&w2)) w2=5;
    /* sliding window width */

    nw = w1*w2;
    win = sf_floatalloc2(w1,w2);

    how = sf_getstring("how"); 
    /* what to compute 
       (fast median, slow median, mean) */
    if (NULL == how) how="fast";

    for (i3=0; i3 < n3; i3++) {

	/* read data plane */
	sf_floatread(data[0],n1*n2,in);
	
	for (i2=0; i2 < n2; i2++) {
	    s2 = SF_MAX(0,SF_MIN(n2-w2,i2-w2/2-1));
	    for (i1=0; i1 < n1; i1++) {
		s1 = SF_MAX(0,SF_MIN(n1-w1,i1-w1/2-1));

		/* copy window */
		for (j2=0; j2 < w2; j2++) {
		    for (j1=0; j1 < w1; j1++) {
			win[j2][j1] = data[s2+j2][s1+j1];
		    }}

		switch (how[0]) {
		    case 'f': /* fast median */
			signal[i2][i1] = 
			    sf_quantile(nw/2,nw,win[0]);
			break;
		    case 's': /* slow median */
			signal[i2][i1] = 
			    slow_median(nw,win[0]);
			break;
		    case 'm': /* mean */
			/* !!! ADD CODE !!! */
			break;
		    default:
			sf_error("Unknown method \"%s\"",how);
			break;
		}
	    }
	}
	
	/* write out */
	sf_floatwrite(signal[0],n1*n2,out);
    }	

    exit(0);
}
