#include <rsf.h>

int main(int argc, char* argv[])
{
  bool adj;
  int nx, ny, nt, N, ix, iy, it, i, j, nw, iw;
  float dx, dy, dt, x0, y0, t0, x, y;
  float *dT, *X, *w;

  sf_file inp, out, gather;

  sf_init(argc, argv);
  inp = sf_input("in");
  out = sf_output("out");
  gather = sf_input("gather");

  /* Adjoint flag */
  if (!sf_getbool("adj",&adj)) adj=false;

  if (adj) {
    /* input dT vector (t^2-t0^2), output w vector (Wx, Wy, Wxy)' */
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&ny)) sf_error("No n2=");
    if (!sf_histint(inp,"n3",&nt)) sf_error("No n3=");
    if (!sf_histfloat(inp,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dy)) sf_error("No d2=");
    if (!sf_histfloat(inp,"d3",&dt)) sf_error("No d3=");
    if (!sf_histfloat(inp,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"o2",&y0)) sf_error("No o2=");
    if (!sf_histfloat(inp,"o3",&t0)) sf_error("No o3=");

    if (!sf_getint("nw",&nw)) nw=3;

    sf_putint(out,"n1",nw);
    sf_putint(out,"n2",nt);
    sf_putint(out,"n3",1);

  } else {
    /* input w vector (Wx, Wy, Wxy), output dT vector (t^2-t0^2) */
    if (!sf_histint(inp,"n1",&nw)) sf_error("No n1=");

    if (!sf_histint(gather,"n1",&nx)) sf_error("No n1=");
    if (!sf_histint(gather,"n2",&ny)) sf_error("No n2=");
    if (!sf_histint(gather,"n3",&nt)) sf_error("No n3=");
    if (!sf_histfloat(gather,"d1",&dx)) sf_error("No d1=");
    if (!sf_histfloat(gather,"d2",&dy)) sf_error("No d2=");
    if (!sf_histfloat(gather,"d3",&dt)) sf_error("No d3=");
    if (!sf_histfloat(gather,"o1",&x0)) sf_error("No o1=");
    if (!sf_histfloat(gather,"o2",&y0)) sf_error("No o2=");
    if (!sf_histfloat(gather,"o3",&t0)) sf_error("No o3=");

  }

  N=nx*ny;
 
  w = sf_floatalloc(nw);
  X = sf_floatalloc(nw);
  dT = sf_floatalloc(N);

  /* Initial Filter */
  if (adj) {
    for(i=0; i < nw; i++) w[i]=0.0;
  } else {
    sf_floatread(w,nw,inp);
  }

  
  /* Loop through t0 coordinate*/
  for (it=0; it < nt; it++) {

    /* Initial Data */
    if (adj) {
      sf_floatread(dT,N,inp);
      for(i=0; i < nw; i++) {
	w[i]=0.0;
      }
    } else {
      for(i=0; i < N; i++) dT[i]=0.0;
    }
    j=0; /*CHANGED FROM JUST ABOVE it LOOP!!!!*/
    /* Loop through x-y coordinates*/
    for (iy=0; iy < ny; iy++){    
      y=2*(y0+iy*dy);

      for (ix=0; ix < nx; ix++){
	x=2*(x0+ix*dx);

	X[0]=x*x;
	X[1]=y*y;
	X[2]=2*x*y;

	/* Loop through slowness parameters*/
	for (iw=0; iw < nw; iw++){

	  /* Solve */
	  if (adj) {
	    w[iw]+=X[iw]*dT[j];
	  } else {
	    dT[j]+=X[iw]*w[iw];
	  }

	}
	/*x-y loops corresponds to loop 0:N-1*/
	j++;
      }
    }
    
    j=0;
    if (adj) sf_floatwrite(w,nw,out);
    if (!adj) sf_floatwrite(dT,N,out);
  }

  exit(0);
}
