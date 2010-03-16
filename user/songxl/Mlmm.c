#include <rsf.h>
#include <math.h>
#include <stdio.h>
int main(int argc, char* argv[]){

    int it, n;
    float *x, x0;
    float t0, h;
    float a;
    sf_file out;
   
     /* initialization */
    sf_init(argc,argv);
    out = sf_output("out");

    /* get trace parameters */

    if (!sf_getfloat("h",&h)) sf_error("Need h=");
    if (!sf_getfloat("seta",&a)) sf_error("Need seta=");
    if (!sf_getfloat("x0",&x0)) x0=1.0;

    t0=0.0;
    n = (int) (1.0/h)+1; 
    sf_warning("n = %d",n);
   /*  change output data dimensions*/
    sf_putint(out,"n1",n+1);
    sf_putfloat(out,"d1",h);
    sf_putfloat(out,"o1",t0);

    x = sf_floatalloc(n+1);

    x[0]=x0;
    x[1]=x[0]+h*x[0];
    x[1]=x[0]+0.5*h*(x[0]+x[1]);
    for (it=2; it< n+1; it++){
        x[it] = -2.0*a*x[it-1]+(1.0+2.0*a)+2.0*h*(1.0+a)*x[it-1];
        }       
 
   sf_floatwrite(x,n+1,out);
   exit(0);
}     
