#include <rsf.h>

#include "stdplot.h"

int main(int argc, char* argv[])
{
    sf_file in;

    sf_init (argc,argv);
    in = sf_input("in");

    vp_erase();

    vp_stdplot_init (0.,2.,0.,1.,
		     true,false,true,true);
    vp_frame_init(in,"blt");
    vp_frame();

    vp_erase();

    vp_stdplot_init (0.,0.036,0.,0.9,
		     true,false,true,false);    
    vp_frame_init(in,"tlb");
    vp_frame();

    exit(0);
}
