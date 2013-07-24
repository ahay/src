/* Create a dipping layer model for HTI testing purposes.  Has fixed velocity structure, but can change dip of layer and degree of anisotropy.*/

/*
  Copyrite (C) 2012 University of Western Australia

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public Licences as published by
  the Free Software Foundation; either version 2 of the License, or 
  (at you option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License 
  along with this program; if not, write to the Free Software Foundation
  Inc., 59 Temle Place, Suite 330, Boston, MA 02111-1307, USA.
*/

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* 
 * input:  Matrix with dims:           uu(nz,nx,ny)
 * output: Stiffness matrix with dims: cc(nz,nx,ny,9)
 */

int main(int argc, char* argv[])
{
	bool verb=true,aniso=false,allaniso=false;

  	sf_file Fu=NULL; /* Input model size */
  	sf_file Fc=NULL; /* Stiffness model */
  	sf_axis a1,a2,a3,a4; /* Cube axes */

  	int   ***mo=NULL;
  	float ***cc11=NULL;
  	float ***cc22=NULL;
  	float ***cc33=NULL;
  	float ***cc44=NULL;
  	float ***cc55=NULL;
  	float ***cc66=NULL;
  	float ***cc12=NULL;
  	float ***cc13=NULL;
  	float ***cc23=NULL;
  	float *vp=NULL,*vs=NULL,*eps=NULL,*del=NULL,*gam=NULL;
  	float *ic11=NULL,*ic33=NULL,*ic44=NULL,*ic55=NULL,*ic13=NULL;
  	float **shit=NULL,**shit2=NULL;

  	int i1,i2,i3;
  	int n1,n2,n3;
  	int lind;

 	int ompchunk;
  	int ompnth=1;
  	int nlayer=6;

  	float rho,ein,din,gin;
  	float d,a,b,x,y,c;
  	float smean;

  	/*----------------------------------------*/
  	/* init RSF */
  	sf_init(argc,argv);
  
  	ompnth=1; /*#omp_init();*/
  	printf("Number of threads: %i",ompnth);


  	if (! sf_getint("ompchunk",&ompchunk)) ompchunk=1; /* set the omp chunk size */
  	if (! sf_getbool("verb",&verb)) verb=true; /* verbose or note (Y/n) */
  	if (! sf_getbool("aniso",&aniso)) aniso=false; /* flag (y/N) for anisotropic layer #2 */
  	if (! sf_getbool("allaniso",&allaniso)) allaniso=false; /* flag (y/N) whether entire model is anisotropic */
  	if (! sf_getfloat("ein",&ein)) ein=.1; /* epsilon anisotropy parameter */
  	if (! sf_getfloat("din",&din)) din=.1; /* delta anisotropy parameter */
  	if (! sf_getfloat("gin",&gin)) gin=.2; /* gamma anisotropy parameter */
  	if (! sf_getfloat("rho",&rho)) rho=2.; /* Background Density model */
	
  	Fu = sf_input ("in") ; /* Input field */
  	Fc = sf_output("out"); /* Output field */

  	a1=sf_iaxa(Fu,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
  	a2=sf_iaxa(Fu,2); sf_setlabel(a2,"a2"); if(verb) sf_raxa(a2);
  	a3=sf_iaxa(Fu,3); sf_setlabel(a3,"a3"); if(verb) sf_raxa(a3);
  	a4=sf_iaxa(Fu,4); sf_setlabel(a4,"a4"); if(verb) sf_raxa(a4);
   
  	sf_oaxa(Fc,a1,1);
  	sf_oaxa(Fc,a2,2);
  	sf_oaxa(Fc,a3,3);
  	sf_oaxa(Fc,a4,4);

  	n1 = sf_n(a1); 
  	n2 = sf_n(a2);
  	n3 = sf_n(a3);
/*  	n4 = sf_n(a4); */

  	/* Model parameter */
  	if (! sf_getfloat("d",&d)) d=n3/6.; /*Parameter in dipping plane: ax+by+cz+d=0 */
  	if (! sf_getfloat("a",&a)) a=1; /*Parameter in dipping plane: ax+by+cz+d=0 */
  	if (! sf_getfloat("b",&b)) b=1; /*Parameter in dipping plane: ax+by+cz+d=0 */
  	c = (a+b)/d*n3;

  	cc11 = sf_floatalloc3(n3,n2,n1);
  	cc22 = sf_floatalloc3(n3,n2,n1);
  	cc33 = sf_floatalloc3(n3,n2,n1);
  	cc44 = sf_floatalloc3(n3,n2,n1);
  	cc55 = sf_floatalloc3(n3,n2,n1);
  	cc66 = sf_floatalloc3(n3,n2,n1);
  	cc12 = sf_floatalloc3(n3,n2,n1);
  	cc13 = sf_floatalloc3(n3,n2,n1);
  	cc23 = sf_floatalloc3(n3,n2,n1);
  	mo = sf_intalloc3(n3,n2,n1);

  	/* Material Properties */
  	vp = sf_floatalloc(nlayer);
  	vs = sf_floatalloc(nlayer);
  	eps= sf_floatalloc(nlayer);
  	del= sf_floatalloc(nlayer);
  	gam= sf_floatalloc(nlayer);

  	/* Stiffness coefficients */
  	ic11 = sf_floatalloc(nlayer);
  	ic33 = sf_floatalloc(nlayer);
  	ic44 = sf_floatalloc(nlayer);
  	ic55 = sf_floatalloc(nlayer);
  	ic13 = sf_floatalloc(nlayer);

  	/* P-wave velocity */
  	vp[0] = 1.5;
  	vp[1] = 2.0;
  	vp[2] = 2.5;
  	vp[3] = 3.0;
  	vp[4] = 3.5;
  	vp[5] = 4.0;
  	for (i1=0; i1<nlayer; i1++){
		vs [i1] = vp[i1]/sqrtf(3.);
		eps[i1] =0.;
		del[i1] =0.;
		gam[i1] =0.;
  	}
  
	/* Construct anisotropic layer */
 	if (aniso) {
		eps[2]=ein;
		del[2]=din;
		gam[2]=gin;
  	}
  
  	if (allaniso){
		for (i1=0; i1 < nlayer; i1++) {
		
			eps[i1]=ein;
			del[i1]=din;
			gam[i1]=gin;  			
			printf("%i,%f,%f,%f\n",i1,eps[i1],del[i1],gam[i1]);
   		}
	}


  	/* Construct material property by layer */
  	for (i1=0; i1 < nlayer; i1++){
		ic33[i1] = rho*vp[i1]*vp[i1];
		ic55[i1] = rho*vs[i1]*vs[i1];
		ic44[i1] = ic55[i1]*(1.+2.*gam[i1]);
		ic11[i1] = ic33[i1]*(1.+2.*eps[i1]);
		ic13[i1] = sqrtf(2.*ic33[i1]*(ic33[i1]-ic55[i1])*del[i1]+(ic33[i1]-ic55[i1])*(ic33[i1]-ic55[i1]))-ic55[i1];
  	}

  	/* Initialize Model */
  	for (i1=0; i1 < n1; i1++) {	
		for (i2=0; i2 < n2; i2++) {
  			for (i3=0; i3 < n3; i3++) {
				mo[i1][i2][i3] = 0;
  			}
		}
  	}

  	/* Second layer*/
  	for (i1=0; i1 < n1; i1++) {	
		for (i2=0; i2 < n2; i2++) {
  			for (i3=floor(n3/6.); i3 < n3; i3++) {
				mo[i1][i2][i3] = 1;
  			}
		}
  	}

  	/* Dipping layer */
  	shit = sf_floatalloc2(n2,n1);
  	shit2= sf_floatalloc2(n2,n1);
	smean = 0.0f;
  	for (i1=0; i1 < n1; i1++) {  x=i1/(n1-1.);	
		for (i2=0; i2 < n2; i2++) {y=i2/(n2-1.);
  			shit[i1][i2] = fminf(n3,fmaxf(1,((d-a*x-b*y)/c)));
  			smean += shit[i1][i2];
		}
  	}
	smean /= (n1*n2);
  
  	/* Scale to desired location */
  	for (i1=0;i1<n1;i1++){
		for (i2=0;i2<n2;i2++){
  			shit2[i1][i2] = (shit[i1][i2] - smean )*400.+2.5*n3/12.;
		}
  	}

  	/* Additional flat layers */
  	for (i1=0; i1 < n1; i1++) {	
		for (i2=0; i2 < n2; i2++) {
  			/* Layer 4 */
  			nlayer=floor(shit2[i1][i2]);
  			for (i3=nlayer; i3 < n3; i3++) {
				mo[i1][i2][i3] = 2;
  			}

  			/* Layer 4 */
  			for (i3=floor(2.*n3/3.); i3 < n3; i3++) {
				mo[i1][i2][i3] = 3;
  			}
  
  			/* Layer 5 */
  			for (i3=floor(9.5*n3/12.); i3 < n3; i3++) {
				mo[i1][i2][i3] = 4;
  			}
  
  			/* Layer 6 */
  			for (i3=floor(11.*n3/12.); i3 < n3; i3++) {
				mo[i1][i2][i3] = 5;
  			}
		}
  	}

  	/* Generate Output model */
  	for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
  			for (i1=0; i1 < n1; i1++) {	
				lind = mo[i1][i2][i3];
				/*printf("%i,%i,%i,%i,%f\n",i1,i2,i3,lind,ic11[lind]);*/
				cc11[i1][i2][i3] = ic11[lind];
				cc22[i1][i2][i3] = ic33[lind];
				cc33[i1][i2][i3] = ic33[lind];
				cc44[i1][i2][i3] = ic44[lind];
				cc55[i1][i2][i3] = ic55[lind];
				cc66[i1][i2][i3] = ic55[lind];
				cc12[i1][i2][i3] = ic13[lind];
				cc13[i1][i2][i3] = ic13[lind];
				cc23[i1][i2][i3] = ic33[lind]-2.*ic44[lind];
  			}
		}
  	}      

  	printf("After final model\n");

  	/* Write output file */
  	sf_floatwrite(cc11[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc22[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc33[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc44[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc55[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc66[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc12[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc13[0][0],n1*n2*n3,Fc);
  	sf_floatwrite(cc23[0][0],n1*n2*n3,Fc);

  	printf("After file writeout\n");

  	/*----------------------------------------*/
  	/* deallocate */
  	free(*cc11);free(cc11);
  	free(*cc22);free(cc22);
  	free(*cc33);free(cc33);
  	free(*cc44);free(cc44);
  	free(*cc55);free(cc55);
  	free(*cc66);free(cc66);
  	free(*cc12);free(cc12);
  	free(*cc13);free(cc13);
  	free(*cc23);free(cc23);
  	free(*mo );free(mo);
  	free(vp);
  	free(vs);
  	free(eps);
  	free(del);
  	free(gam);
  	free(ic11);
  	free(ic33);
  	free(ic44);
  	free(ic55);
  	free(ic13);
  	free(*shit);free(shit);
  	free(*shit2);free(shit2);
  	/*----------------------------------------*/

  	exit(0);
}
