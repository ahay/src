/* Apply running mean or median filter */

#include <rsf.h>

static float median(int n, float *a)
/* find median by the fast quantile, changes a */
{
    float *i, *j, ak, *low, *hi, buf, *k;

    low=a;
    hi=a+n-1;
    k=a+n/2; 
    while (low<hi) {
	ak = *k;
	i = low; j = hi;
	do {
	    while (*i < ak) i++;     
	    while (*j > ak) j--;     
	    if (i<=j) {
		buf = *i;
		*i++ = *j;
		*j-- = buf;
	    }
	} while (i<=j);
	if (j<k) low = i; 
	if (k<i) hi = j;
    }
    return (*k);
}

int main(int argc, char* argv[]) 
{
    int w1, w2, nw, shift1,shift2, j1,j2, i1,i2,i3, n1,n2,n3;
    int what;
    float **data, **signal, **win;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* get data dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    /* input and output */
    data = sf_floatalloc2(n1,n2);
    signal = sf_floatalloc2(n1,n2);

    if (!sf_getint("w1",&w1)) w1=5;
    if (!sf_getint("w2",&w2)) w2=5;
    /* sliding window width */

    nw = w1*w2;
    win = sf_floatalloc2(w1,w2);

    if (!sf_getint("what",&what)) what=0;
    /* [0,1,2] 0: median, 1: mean, 2: LUM */

    for (i3=0; i3 < n3; i3++) {

	/* read data plane */
	sf_floatread(data[0],n1*n2,in);
	
	for (i2=0; i2 < n2; i2++) {
	    shift2 = SF_MAX(0,SF_MIN(n2-w2,i2-w2/2-1));
	    for (i1=0; i1 < n1; i1++) {
		shift1 = SF_MAX(0,SF_MIN(n1-w1,i1-w1/2-1));

		/* copy window */
		for (j2=0; j2 < w2; j2++) {
		for (j1=0; j1 < w1; j1++) {
		    win[j2][j1] = data[shift2+j2][shift1+j1];
		}}

		switch (what) {
		    case 0: /* median */
			signal[i2][i1] = median(nw,win[0]);
			break;
		    case 1: /* mean */
		    case 2: /* LUM  */
		    default:
			/* ADD CODE */
			break;
		}
	    }
	}
	
	/* write out */
	sf_floatwrite(signal[0],n1*n2,out);
    }	

    exit(0);
}
