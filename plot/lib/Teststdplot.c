#include <rsf.h>

#include "stdplot.h"

int main(int argc, char* argv[])
{
    sf_file in;

    sf_init (argc,argv);
    in = sf_input("in");

    vp_stdplot_init (0.,1.,0.,2.);
    vp_frame_init(in,7);
    vp_frame();

    exit(0);
}
