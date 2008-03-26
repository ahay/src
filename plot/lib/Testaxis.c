#include <stdio.h>

#include "axis.h"

int main (void)
{
    int n, maxstrlen;
    float onum, dnum;

    n = vp_optimal_scale(10,true,true,"%1.5g",0.,100.,&onum,&dnum,&maxstrlen);
    printf("min/max: 0./100.: %d steps of %g from %g to %g\n",
	   n-1,dnum,onum,onum+(n-1)*dnum);

    n = vp_optimal_scale(10,true,true,"%1.5g",0.95,2.01,&onum,&dnum,&maxstrlen);
    printf("min/max: 0.95/2.01: %d steps of %g from %g to %g\n",
	   n-1,dnum,onum,onum+(n-1)*dnum);
    
    n = vp_optimal_scale(10,true,true,"%1.5g",
			 0.,20.0008,&onum,&dnum,&maxstrlen);
    printf("min/max: 0./20.0008: %d steps of %g from %g to %g\n",
	   n-1,dnum,onum,onum+(n-1)*dnum);

    n = vp_optimal_scale(10,true,true,"%1.5g",
			 -2.89735e-17,0.64849,&onum,&dnum,&maxstrlen);
    printf("min/max: -2.89735e-17/0.64849: %d steps of %g from %g to %g\n",
	   n-1,dnum,onum,onum+(n-1)*dnum);

    n = vp_optimal_scale(10,true,true,"%1.5g",
			 -0.680925,-0.651386,&onum,&dnum,&maxstrlen);
    printf("min/max: -0.680925/-0.651386: %d steps of %g from %g to %g\n",
	   n-1,dnum,onum,onum+(n-1)*dnum);

    return 0;
}
