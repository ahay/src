#include <rsf.h>

int main(int argc, char* argv[])
{
  bool adj;
  int nx, nt, ix, it, i, nm, im;
  float dx, dt, x0, t0, x;
  float *dT=NULL, *X=NULL, *m=NULL;

  sf_file inp=NULL, out=NULL, gather=NULL;

  sf_init(argc, argv);
  inp = sf_input("in");
  out = sf_output("out");
  gather = sf_input("gather");

  /* Adjoint flag */
  if (!sf_getbool("adj",&adj)) adj=false;

  if (adj) {
    /* input dT vector x^2/ (t^2-t0^2), output m vector (v^2_nmo, a)' */
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nt)) sf_error("No n2=");
    if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dt)) sf_error("No d2=");
    if (!sf_histfloat(inp,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"o2",&t0)) sf_error("No o2=");

    if (!sf_getint("nm",&nm)) nm=2;

    sf_putint(out,"n1",nm);
    sf_putint(out,"n2",nt);
    sf_putint(out,"n3",1);

  } else {
    /* input m vector (v^2_nmo, a), output dT vector x^2/(t^2-t0^2) */
    if (!sf_histint(inp,"n1",&nm)) sf_error("No n1=");

    if (!sf_histint(gather,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(gather,"n2",&nt)) sf_error("No n2=");
    if (!sf_histfloat(gather,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(gather,"d2",&dt)) sf_error("No d2=");
    if (!sf_histfloat(gather,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(gather,"o2",&t0)) sf_error("No o2=");

  }

  m = sf_floatalloc(nm);
  X = sf_floatalloc(nm);
  dT = sf_floatalloc(nx);

  /* Initial Filter */
  if (adj) {
    for(i=0; i < nm; i++) m[i]=0.0;
  } else {
    sf_floatread(m,nm,inp);
  }

  /* Loop through t0 coordinate*/
  for (it=0; it < nt; it++) {

    /* Initial Data */
    if (adj) {
      sf_floatread(dT,nx,inp);
      for(i=0; i < nm; i++) {
	m[i]=0.0;
      }
    } else {
      for(i=0; i < nx; i++) dT[i]=0.0;
    }

    /* Loop through x-coordinates*/

    for (ix=0; ix < nx; ix++){
	x=(x0+ix*dx);

	X[0]=1.0;
	X[1]=x*x;

	/* Loop through model parameters*/
	for (im=0; im < nm; im++){

	  /* Solve */
	  if (adj) {
	    m[im]+=X[im]*dT[ix];
	  } else {
	    dT[ix]+=X[im]*m[im];
	  }
      }
    }

    if (adj) sf_floatwrite(m,nm,out);
    if (!adj) sf_floatwrite(dT,nx,out);
  }

  exit(0);
}
