#include <rsf.h>

#include "random.h"


int main(int argc, char* argv[])
{
    int n1, n2, i2, i, j;
    float *data, sigma, rand, thresh, thresh2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getfloat("sigma",&sigma)) sigma=1.;
    /* noise magnitude */

    if (!sf_getfloat("thresh",&thresh)) thresh=0.93;
    /* noise threshold */

    if (!sf_getfloat("thresh2",&thresh2)) thresh2=0.4;
    /* noise threshold */

    data = sf_floatalloc(n1);

    random_init (121094);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	i=0;
	while (i < n1) {
	    rand = random0();
	    if (rand < thresh) { /* no noise */
		i++;
	    } else { 
		for (j=i; j < n1; j++) {
		    rand = random0();
		    data[j] += sigma*rand;
	
		    rand = random0();
		    if (rand < thresh2) break;
		}
		i = j + 1;
	    }
	}
	
	sf_floatwrite(data,n1,out);
    }

/*
		call getnoise(sigma, noise)
		data = data + noise
*/

    sf_close();
    exit(0);
}
