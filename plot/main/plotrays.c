#include <rsf.h>
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    int n1, n2, ir, jr;
    float o1, d1, o2, d2, **traj;
    int nsr, it, nt;

    sf_init (argc,argv);
    vp_init ();

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    if (!sf_getint("n2",&n2)) sf_error("Need n2=");
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    if (!sf_getfloat("d2",&d2)) sf_error("Need d2=");
    if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");
    if (!sf_getfloat("o2",&o2)) sf_error("Need o2=");
    if (!sf_getint("nt",&nt)) nt=n1*n2;
    if (!sf_getint("jr",&jr)) jr=1;

    traj = sf_floatalloc2 (2,nt);

    /* transp and yreverse */
    vp_stdplot_init (o1,o1+(n1-1)*d1,o2,o2+(n2-1)*d2,
		     true,false,true,false);
    vp_plot_init(1);
    vp_plot_set(0);
 
    fread(&nsr,sizeof(int),1,stdin);

    for (ir=0; ir < nsr; ir++) {
	fread(&it,sizeof(int),1,stdin);
	fread(traj[0],sizeof(float),(it+1)*2,stdin);
	if (ir>0 && ir%jr) continue;  
	vp_umove(traj[it][1],traj[it][0]);
	while (--it >= 0) {
	    vp_udraw(traj[it][1],traj[it][0]);
	}
    }

    vp_frame();

    exit(0);
}
