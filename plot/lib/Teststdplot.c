#include <rsf.h>

#include "stdplot.h"

int main(int argc, char* argv[])
{
    sf_file in;

    sf_init (argc,argv);
    in = sf_input("in");

    vp_erase();

    vp_stdplot_init (0.,1.,2.,0.);
    vp_frame_init(in, true);
    vp_frame();

    vp_erase();

    vp_stdplot_init (0.,0.9,0.036,0.);
    vp_frame_init(in, true);
    vp_frame();

    exit(0);
}
