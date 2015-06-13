#include <rsf.h>

#include "cosft3.h"

int main(int argc, char*argv[])
{
    bool inv;
    int n1, n2, n3;
    float ***dat; 
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&n3)) sf_error("No n3=");

    if (!sf_getbool("inv",&inv)) inv=false;

    dat = sf_floatalloc3(n1,n2,n3);

    sf_floatread(dat[0][0],n1*n2*n3,inp);
    cosft3(inv,n1,n2,n3,dat);
    sf_floatwrite(dat[0][0],n1*n2*n3,out);

    exit(0);
}
