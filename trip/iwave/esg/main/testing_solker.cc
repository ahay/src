
#include <stdio.h>
#include <stdlib.h>
// #include <ansol_esgsteps.h>

int main(void){

	int indexes[3];
		indexes[0] = 1;
		indexes[1] = 0;
		indexes[2] = 1;
	int ndim = 2;
	float x[2]; 
		x[0] = 1.1;  
		x[1] = 2.22;
	float xs[2];
		xs[0] = -4.0; 
		xs[1] = 0.0;
	float t = 3.11;
	float alpha = 2;
	float beta = 3;
	float fpeak = 0.05;

	float out_primexk, out_primet;

// 	out_primexk = uij_primexk( indexes,
// 		     		   ndim,
// 		    		   x,
// 		     		   xs,
// 				   t,
// 				   alpha,
// 				   beta,
// 			  	   fpeak,
// 			 	   stderr );
// 
// 	out_primet = uij_primet( indexes,
// 		     		   ndim,
// 		    		   x,
// 		     		   xs,
// 				   t,
// 				   alpha,
// 				   beta,
// 			  	   fpeak,
// 			 	   stderr );
// 
// 	fprintf(stderr,"out_primexk = %g\n",out_primexk);
// 	fprintf(stderr,"out_primet  = %g\n",out_primet);

	return 0;
}