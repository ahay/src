#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int n, n1, n2, i1, i2;
    float d, g, s, v0, v, x, z, a;
    float **vel, **time;
    sf_file out;

    /* SEPlib initialization */
    sf_init (argc,argv);
    out = sf_output("out");

    if (!sf_getint("n",&n)) sf_error("Need n=");
    if (!sf_getfloat("g",&g)) g = 1.;
    if (!sf_getfloat("v0",&v0)) v0 = 0.5;
    if (!sf_getfloat("s",&s)) s = 0.5;
    
    d = 0.5/n;
    n1 = n+1;
    n2 = 2*n+1;

    a = 0.5*g*g/v0;

    sf_putint(out,"n1",n1); sf_putfloat(out,"d1",d); sf_putfloat(out,"o1",0.);
    sf_putint(out,"n2",n2); sf_putfloat(out,"d2",d); sf_putfloat(out,"o2",0.);
    sf_putint(out,"n3",2);
    sf_setformat(out,"native_float");

    vel = sf_floatalloc2 (n1,n2);
    time = sf_floatalloc2 (n1,n2);

    for (i2 = 0; i2 < n2; i2++) {
      x = i2*d - s;
      x = x*x;
      for (i1 = 0; i1 < n1; i1++) {
	z = i1*d;
	v = v0 + g*z;
	z = 1. + a*(z*z+x)/v;
	vel[i2][i1] = v;
	time[i2][i1] = log(z + sqrt(z*z-1.))/g;
      }
    }

    sf_write(vel[0],sizeof(float),n1*n2,out);
    sf_write(time[0],sizeof(float),n1*n2,out);

    exit (0);
}

  

 
