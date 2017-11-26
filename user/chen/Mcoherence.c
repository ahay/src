/* 3D Coherence Cube, C1, C2, C3 in one */

/*
  Copyright (C) 2013 Zhonghuan Chen, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranyw of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

float coh1(float **u, int n1, int n2);
float coh2(float **u, int n1, int n2);
float coh3(float **u, int n1, int n2);

int main(int argc, char* argv[])
{
    sf_file in, out, idip, xdip;
    int n1, n2, n3, np1, np2, ntw, nxw, nyw, m1, m2;
    int i1, i2, i3, j1, j2, j3, k1, k2, ix, iy;
    float ***u1, ***u2, **u3, *p, it;
    float coh, max, ***v1=NULL, ***v2=NULL, d1, d2;
    bool twod, verb;
    char *mode;
    float (*coh_fun)(float **, int , int);
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


    switch(mode[1])
    {
	case '2': coh_fun=coh2; break;
	case '3': coh_fun=coh3; break;
	default: coh_fun=coh1; nxw=1; nyw=1;
    }
    if(n3==1)  twod=true;
    if(twod) {nyw=0; np2=1; op2=0.0;}

    u1 = sf_floatalloc3(n1, n2, n3);
    u2 = sf_floatalloc3(n1, n2, n3);
    if(idip)v1 = sf_floatalloc3(n1, n2, n3);
    if(xdip)v2 = sf_floatalloc3(n1, n2, n3);

    m1 = 2*ntw+1;
    m2 = (2*nxw+1)*(2*nyw+1);
    u3 = sf_floatalloc2(m1, m2);

    sf_floatread(u1[0][0], n1*n2*n3, in);
    for(i3=0; i3<n3; i3++)
    {
	for(i2=0; i2<n2; i2++)
	    for(i1=0; i1<n1; i1++)
	    {
		max = 0.0;
		d1 = 0.0;
		d2 = 0.0;
		for(k2=0; k2<np2; k2++)
		    for(k1=0; k1<np1; k1++)
		    {
			p = u3[0];
			p1 = op1+dp1*k1;
			p2 = op2+dp2*k2;
			for(j3=-nyw; j3<=nyw; j3++)
			    for(j2=-nxw; j2<=nxw; j2++)
				for(j1=-ntw; j1<=ntw; j1++, p++)
				{
				    iy = i3+j3;
				    ix = i2+j2;
				    it = i1+j1+j3*p2+j2*p1;
				    if( iy<0 || iy>=n3 || ix<0 || ix>=n2 || it<0 || it>=n1-1)
					*p = 0.0;
				    else *p = (((int)it)+1-it)*u1[iy][ix][(int)it]+
					     (it-(int)it)*u1[iy][ix][(int)it+1];
				}
			coh = coh_fun(u3, m1, m2);
			if(coh>max) 
			{
			    max=coh;
			    d1 = p1;
			    d2 = p2;
			}
		    }
		u2[i3][i2][i1] = max;
		if(idip) v1[i3][i2][i1] = d1;
		if(xdip) v2[i3][i2][i1] = d2;
	    }
	if(verb) sf_warning("%d of %d;", i3, n3);
    }
    sf_floatwrite(u2[0][0], n1*n2*n3, out);
    if(idip)sf_floatwrite(v1[0][0], n1*n2*n3, idip);
    if(xdip)sf_floatwrite(v2[0][0], n1*n2*n3, xdip);

    free(**u1);
    free(*u1);
    free(u1);
    free(**u2);
    free(*u2);
    free(u2);
    free(*u3);
    free(u3);
    if(idip)
    {
	free(**v1);
	free(*v1);
	free(v1);
    }
    if(xdip)
    {
	free(**v2);
	free(*v2);
	free(v2);
    }
    return 0;
}

float coh1(float **u, int n1, int n2)
{
    float *p1, *p2, *p3, c11, c22, c33, c12, c13;
    int i1;

    if(n2==3) // 2D
    {
	p1 = u[1];
	p2 = u[0];
	c11 = 0.0;
	c22 = 0.0;
	c12 = 0.0;
	for(i1=0; i1<n1; i1++)
	{
	    c12 += p1[i1]*p2[i1];
	    c11 += p1[i1]*p1[i1];
	    c22 += p2[i1]*p2[i1];
	}
	c11 = sqrt(c11*c22);
	if(c11 == 0.0) return 1.0;
	else return (c12/c11);
    }else if(n2 == 9) // 3D
    {
	p1 = u[4];
	p2 = u[5];
	p3 = u[1];
	c11 = 0.0;
	c22 = 0.0;
	c33 = 0.0;
	c12 = 0.0;
	c13 = 0.0;
	for(i1=0; i1<n1; i1++)
	{
	    c11 += p1[i1]*p1[i1];
	    c22 += p2[i1]*p2[i1];
	    c33 += p3[i1]*p3[i1];
	    c12 += p1[i1]*p2[i1];
	    c13 += p1[i1]*p3[i1];
	}
	c11 *= sqrt(c22*c33);
	if(c11 == 0.0) return 1.0;
	return sqrt(c12*c13/c11);
    }else return 1.0;
}


float coh2(float **u, int n1, int n2)
{
    float d1, d2, coh;
    int i1, i2;

    coh = 0.0;
    d2 = 0.0;
    for(i1=0; i1<n1; i1++)
    {
	d1 = 0.0;
	for(i2=0; i2<n2; i2++)
	{
	    d1 += u[i2][i1];
	    d2 += u[i2][i1] * u[i2][i1];
	}
	coh += d1*d1;
    }
    return (coh/d2/n2);
}


float coh3(float **u, int n1, int n2)
{
    float d1, d2, *work, *s; /* **c, alpha=1.0, beta=0.0; */
    int i1, ns, lwork, info;

    lwork = n1*n2;
    work = sf_floatalloc(lwork);
    ns = n1<n2?n1:n2;
    s = sf_floatalloc(lwork);

// THIS WAY IT CAN NOT CONVERGED !!!!
// =============================
//	c = sf_floatalloc2(n2, n2);
//	ssyrk_("U", "T", &n2, &n1, &alpha, u, &n1, &beta, c[0], &n2);
//	ssyev_("N", "U", &n2, c[0], &n2, s, work, &lwork, &info);
//	if(info!=0) sf_warning("ssss %d", info);
//	free(*c);
//	free(c);

    sgesvd_("N", "N", &n1, &n2, *u, &n1, s, 
	    NULL, &n1, NULL, &n2, work, &lwork, &info);
    free(work);

    d1 = s[0]*s[0]; d2 =0.0;
    for(i1=0; i1<ns; i1++) d2 += s[i1]*s[i1];

    free(s);
    return (d1/d2);
}



