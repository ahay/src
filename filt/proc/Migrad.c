/* Gradient on the first axis. */

#include <rsf.h>
#include "igrad1.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2;
    float *pp, *qq;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    pp = sf_floatalloc(n1);
    qq = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_read(pp,sizeof(float),n1,in);
	igrad1_lop (false,false,n1,n1,pp,qq);
	sf_write(qq,sizeof(float),n1,out);
    }

    exit(0);
}


