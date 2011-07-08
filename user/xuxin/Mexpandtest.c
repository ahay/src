/* testing wave.c */
/* cut not working */
#include <rsf.h>
#include "wave.h"

int main(int argc, char *argv[])
{
    bool fsrf;
    sf_file Fi,Fo,Fc;
    sf_axis ax,az,at;
    wave2dp W;
    float **in=NULL,**ou=NULL,**cu=NULL;
    int nx,nbx,mx,nz,nbz,mz;

    sf_init(argc,argv);
    Fi = sf_input("in");
    Fo = sf_output("out");
    Fc = sf_output("cut");

    if (!sf_histint(Fi,"n1",&nx))   sf_error("No n1= in input");
    if (!sf_histint(Fi,"n2",&nz))   sf_error("No n2= in input");
    if (!sf_getint("nb1",&nbx))     nbx=5;
    if (!sf_getint("nb2",&nbz))     nbz=5;
    if (!sf_getbool("free",&fsrf))  fsrf=false;
    /* if y, cut; if n, expand */

    ax = sf_iaxa(Fi,1); 
    az = sf_iaxa(Fi,2);
    at = sf_maxa(1,0.0,1.0);

    W = wave2d_init(ax,az,at,nbx,nbz,fsrf);
    mx = sf_n(W->px);
    mz = sf_n(W->pz);

    in = sf_floatalloc2(nx,nz);
    ou = sf_floatalloc2(mx,mz);
    cu = sf_floatalloc2(nx,nz);

    sf_floatread(in[0],nx*nz,Fi);

    sf_warning("n1=%d n2=%d -> n1=%d n2=%d",nx,nz,mx,mz);

    /* expand */
    expand2(in,ou,W);
    sf_oaxa(Fo,W->px,1);
    sf_oaxa(Fo,W->pz,2);
    sf_floatwrite(ou[0],mx*mz,Fo);

    /* cut */
    cut2(ou,cu,W);
    sf_oaxa(Fc,ax,1);
    sf_oaxa(Fc,az,2);
    sf_floatwrite(cu[0],nx*nz,Fc);

    exit(0);
}
