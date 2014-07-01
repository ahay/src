/* 3D coherence computation 	
*/

/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
   
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
   
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  References: 
    1) Marfurt, Kurt J., et al. "Coherency calculations in the presence 
	of structural dip." Geophysics 64.1 (1999): 104-111.
    2) Wang, Xiaokai, Jinghuai Gao, Wenchao Chen, and Yongzhong Song. 
	"An efficient implementation of eigenstructure-based coherence 
	algorithm using recursion strategies and the power method." Journal
	 of Applied Geophysics 82 (2012): 11-18.
*/

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif


float coh1(float **cxy, int J)
/*< 1st generation coherence >*/
{
	float t=cxy[0][0]*sqrtf(cxy[1][1]*cxy[2][2]);
	return sqrtf(cxy[0][1]*cxy[0][2]/(t+SF_EPS));
}

float coh2(float **cxy, int J)
/*< 2nd generation coherence >*/
{
	int i,j;
	float s1,s2;
	s1=s2=0;
	for(i=0; i<J; i++)
	{
		for(j=0; j<J; j++) s1+=cxy[i][j];
		s2+=cxy[i][i];
	}
	return s1/(s2+SF_EPS)/J;
}


static	float *u, *v;
float coh3(float **cxy, int J)
/*< find the maximum eigenvalue using power method
NB: sum of eigenvalues=trace of cxy >*/
{
	int i, j, k, maxidx, niter=30;
	float s, t, m1, m;

	s=m1=0;
	for(i=0; i<J; i++) {
		s+=cxy[i][i];/* trace of matrix */
		u[i]=1.0;/* initialize u */
	}
	
	for(k=0; k<niter; k++){
		for(i=0; i<J; i++){
			t=0.0;
			for(j=0; j<J; j++)
				t+=cxy[i][j]*u[j];/* v=A u*/
			v[i]=t;			
		}

		m=fabsf(v[0]); maxidx=0;
		for(i=0; i<J; i++) { 
			maxidx=(m>fabsf(v[i]))?maxidx:i;
			m=(m>fabsf(v[i]))?m:fabsf(v[i]); 
		}
		m=fabsf(v[maxidx]);
		for(i=0; i<J; i++) u[i]=v[i]/m;

		if(fabsf(m-m1)<1.e-6) break;
		m1=m;	
	}
	return m/(s+SF_EPS);
}

int main(int argc, char *argv[])
{
    	sf_file in, out;
    	int n1, n2, n3, ntw, nxw, nyw;
    	int i1, i2, i3, j1, j2, j3, k2, k3, px, py, J;
	float ***u1, ***u2, **cxy;
	char *mode;
	float (*cohn)(float **, int);

    	sf_init(argc, argv);
    	in=sf_input("in");	/* 3D seismic data volume */
   	out=sf_output("out");	/* 3D coherence volume */

    	if (!sf_histint(in,"n1",&n1)) 	sf_error("No n1= in input");
    	if (!sf_histint(in,"n2",&n2)) 	sf_error("No n2= in input");
    	if (!sf_histint(in,"n3",&n3)) 	n3=1;	/* default: n3=1 if 2D */

    	if (!sf_getint("ntw",&ntw)) 	ntw=10; /* radius of the window in t */
    	if (!sf_getint("nxw",&nxw)) 	nxw=1; /* radius of the window in x */
    	if (!sf_getint("nyw",&nyw)) 	nyw=1; /* radius of the window in y */
	if(!(mode=sf_getstring("mode"))) mode = "c3"; /* coherence type: c1,c2,c3 */
	J=(2*nxw+1)*(2*nyw+1);
	switch(mode[1])
	{
		case '3': cohn=coh3; 
			u=(float*)malloc(J*sizeof(float));
			v=(float*)malloc(J*sizeof(float));
			memset(u, 0, J*sizeof(float));
			memset(v, 0, J*sizeof(float));
			break;
		case '2': cohn=coh2; 
			break;
		default: cohn=coh1;
	}

	cxy=sf_floatalloc2(J,J);
	u1 = sf_floatalloc3(n1, n2, n3);
	u2 = sf_floatalloc3(n1, n2, n3);
	sf_floatread(u1[0][0], n1*n2*n3, in);	

	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
	    memset(cxy[0],0,J*J*sizeof(float));
	    for(j1=-ntw; j1<=ntw; j1++) {/* construct the covariance matrix */
		if (i1+j1>=0 && i1+j1<n1)
		{
			for(j3=-nyw; j3<=nyw; j3++) for(j2=-nxw; j2<=nxw; j2++)
			for(k3=-nyw; k3<=nyw; k3++) for(k2=-nxw; k2<=nxw; k2++)
			{
				px=j2+nxw+(2*nxw+1)*(j3+nyw);
				py=k2+nxw+(2*nxw+1)*(k3+nyw);
				if ( 	(i2+j2>=0 && i2+j2<n2) && 
					(i3+j3>=0 && i3+j3<n3) &&
					(i2+k2>=0 && i2+k2<n2) && 
					(i3+k3>=0 && i3+k3<n3) )
				cxy[py][px]+=u1[i3+j3][i2+j2][i1+j1]*u1[i3+k3][i2+k2][i1+j1];
			}
		}
	    }

	    u2[i3][i2][i1]=cohn(cxy, J);
	}
	sf_floatwrite(u2[0][0], n1*n2*n3, out);

	free(**u1); free(*u1); free(u1);
	free(**u2); free(*u2); free(u2);
	free(*cxy); free(cxy);
	if(strcmp(mode,"c3")==0) { free(u); free(v);}

    	exit(0);
}

