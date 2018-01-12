/* Coherence calculations in the presence of structural dip
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
  NB: modified according to ../chen/Mcoherence.c by Zhonghuan Chen.
*/

#include <rsf.h>

float coh1(float **cpq, int J)
/*< 1st generation coherence >*/
{
	float tmp;
	if(J==3){// 2D
		return cpq[0][1]/sqrtf(cpq[0][0]*cpq[1][1]+SF_EPS);
	}else if(J == 9) {// 3D
		tmp=cpq[0][0]*sqrtf(cpq[1][1]*cpq[2][2]);
		return sqrtf(cpq[0][1]*cpq[0][2]/(tmp+SF_EPS));
	}else return 1.0;
}

float coh2(float **cpq, int J)
/*< 2nd generation coherence >*/
{
	int i,j;
	float s1,s2;
	s1=s2=0;
	for(i=0; i<J; i++)
	{
		for(j=0; j<J; j++) s1+=cpq[i][j];
		s2+=cpq[i][i];
	}
	return s1/(s2+SF_EPS)/J;
}

static	float *u, *v;
float coh3(float **cpq, int J)
/*< find the maximum eigenvalue using power method >*/
{
	int i, j, k, maxidx, niter=30;
	float s, t, m1, m;

	s=m1=0;
	for(i=0; i<J; i++) {
		s+=cpq[i][i];/* trace of matrix */
		u[i]=1.0;/* initialize u */
	}
	
	for(k=0; k<niter; k++){
		for(i=0; i<J; i++){
			t=0.0;
			for(j=0; j<J; j++) t+=cpq[i][j]*u[j];/* v=A u*/
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

int main(int argc, char* argv[])
{
	sf_file in, out, idip, xdip;
	int n1, n2, n3, np1, np2, ntw, nxw, nyw, m, J;
	int i1, i2, i3, j1, j2, j3, k1, k2, ix, ia, iy, i, j, k;
	float ***u1, ***u2, **u3, ***v1, ***v2, **cpq;
	float mmax, it, s, d1, d2;
	bool twod, verb;
	char *mode;
	float (*cohn)(float **cpq, int J);
	float op1, op2, dp1, dp2, p1, p2;

	sf_init(argc, argv);

	in  = sf_input("in");	/* 3D data set */
	out = sf_output("out");	/* 3D coherence cube */

	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histint(in,"n3",&n3)) n3=1; 

	if (!sf_getint("ntw", &ntw) ) ntw=5;
	/* half window size for coherence */
	if (!sf_getint("nxw", &nxw) ) nxw=5;
	/* half window size for coherence */
	if (!sf_getint("nyw", &nyw) ) nyw=5;
	/* half window size for coherence */
	if (!sf_getbool("twod",&twod) ) twod=true;
	/* y: only twod coherence */
	if (!sf_getbool("verb",&verb) ) verb=true;
	/* verbosity */
	if(sf_getstring("idip")!=NULL) idip  = sf_output("idip");
	/* inline dip */
	else idip = NULL;
	if(sf_getstring("xdip")!=NULL) xdip  = sf_output("xdip");
	/* crossline dip */
	else xdip = NULL;
	if((mode=sf_getstring("mode"))==NULL) mode = "c1";
	/* coherence mode: c1, c2, c3 */
	if (!sf_getfloat("op1", &op1) ) op1=-2.0;
	if (!sf_getfloat("dp1", &dp1) ) dp1=0.5;
	if (!sf_getint("np1", &np1) ) np1=9;
	/* inline slope */
	if (!sf_getfloat("op2", &op2) ) op2=-2.0;
	if (!sf_getfloat("dp2", &dp2) ) dp2=0.5;
	if (!sf_getint("np2", &np2) ) np2=9;
	/* xline slope */

	m = 2*ntw+1;
	J = (2*nxw+1)*(2*nyw+1);
	switch(mode[1])
	{
		case '3': cohn=coh3; 
			u=(float*)malloc(J*sizeof(float));
			v=(float*)malloc(J*sizeof(float));
			memset(u, 0, J*sizeof(float));
			memset(v, 0, J*sizeof(float));
			break;
		case '2': cohn=coh2; break;
		default: cohn=coh1; nxw=1; nyw=1;
	}
	if(n3==1)  twod=true;
	if(twod) {nyw=0; np2=1; op2=0.0;}

	u1 = sf_floatalloc3(n1, n2, n3);
	u2 = sf_floatalloc3(n1, n2, n3);
	if(idip) v1=sf_floatalloc3(n1, n2, n3);
	if(xdip) v2=sf_floatalloc3(n1, n2, n3);
	u3=sf_floatalloc2(m, J);
	cpq=sf_floatalloc2(J,J);

	sf_floatread(u1[0][0], n1*n2*n3, in);
	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		mmax = 0.0;
		d1 = 0.0;
		d2 = 0.0;
		for(k2=0; k2<np2; k2++)
		for(k1=0; k1<np1; k1++)
		{
			p1 = op1+dp1*k1;
			p2 = op2+dp2*k2;
			memset(u3[0],0,m*J*sizeof(float));
			for(j3=-nyw; j3<=nyw; j3++)
			for(j2=-nxw; j2<=nxw; j2++)
			for(j1=-ntw; j1<=ntw; j1++)
			{
				iy = i3+j3;
				ix = i2+j2;
				it = i1+j1+j3*p2+j2*p1;
				if( iy<0 || iy>=n3 || ix<0 || ix>=n2 || it<0 || it>=n1-1){
					u3[j2+nxw+(2*nxw+1)*(j3+nyw)][j1+ntw] = 0.0;
				}else {//linear interpolation
					ia=(int)it;
					u3[j2+nxw+(2*nxw+1)*(j3+nyw)][j1+ntw]=(ia+1.0-it)*u1[iy][ix][ia]+(it-ia)*u1[iy][ix][ia+1];
				}
			}
			/* construct the covariance matrix */
			for(i=0; i<J; i++)
			for(j=0; j<J; j++)
			{
				s=0;
				for(k=0; k<m; k++) s+=u3[i][k]*u3[j][k];
				cpq[i][j]=s;
			}
			s=cohn(cpq, J);
			if(s>mmax) 
			{
				mmax=s;
				d1 = p1;
				d2 = p2;
			}
		}
		u2[i3][i2][i1] = mmax;
		if(idip) v1[i3][i2][i1] = d1;
		if(xdip) v2[i3][i2][i1] = d2;
	}
	sf_floatwrite(u2[0][0], n1*n2*n3, out);
	if(idip) sf_floatwrite(v1[0][0], n1*n2*n3, idip);
	if(xdip) sf_floatwrite(v2[0][0], n1*n2*n3, xdip);

	free(**u1); free(*u1); free(u1);
	free(**u2); free(*u2); free(u2);
	free(*u3); free(u3);
	free(*cpq); free(cpq);
	if(idip){ free(**v1); free(*v1); free(v1);}
	if(xdip){ free(**v2); free(*v2); free(v2);}
	if(strcmp(mode,"c3")==0) { free(u); free(v);}

	exit(0);
}

