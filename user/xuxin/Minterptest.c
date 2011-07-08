/* testing 1D interpolation */

#include <rsf.h>
#include "interp.h"

#define NX 4
#define NY 3

int main(int argc, char* argv[])
{
    float x[]={-1,0,1,2};
    float fx[]={-1,1,3,5};
    float y[]={-2,0.5,2.5};
    float *fy;
    lint1dp L;
    int i;

    fy = sf_floatalloc(NY);
    L = lint1d_init(x,NX,y,NY);
    lint1d_extract(fx,fy,L);
    
    for (i=0; i < NY; i++)
	printf("fy[%d] = %g\n",i,fy[i]);
  
}
