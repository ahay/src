#include <rsf.h>

#include "stdplot.h"
#include "vplot.h"

int main(int argc, char* argv[])
{
    sf_file in;

    sf_init (argc,argv);
    in = sf_input("in");

    vp_erase();

    vp_stdplot_init (0.,2.,0.,1.,
		     true,false,true,true);
    vp_frame_init(in,"blt");
    vp_framenum(0.0001);
    vp_frame();

    vp_erase();

    vp_stdplot_init (0.,0.036,0.,0.9,
		     true,false,true,false);    
    vp_frame_init(in,"tlb");
    vp_framenum(15.);
    vp_frame();


    vp_barframe_init (0.,1.);
    vp_barframe ();

    exit(0);
}
