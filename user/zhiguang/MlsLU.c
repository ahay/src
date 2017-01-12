/* Local similarity filter (direct solving) */
/*
  Copyright (C) 2016 University of Texas at Austin
   
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
*/
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX(a,b) ((a)>(b)? a : b)
#define MIN(a,b) ((a)<(b)? a : b)
#define SWAP(a,b) {dum=(a); (a)=(b); (b)=dum;}
#define TINY 1.e-18

void bandAx(float **a, int n, int m1, int m2, float *x, float *b)
/*< band matrix Ax=b >*/
{
	int i, j, k, tmploop;

	for (i=1; i<=n; i++){
		k=i-m1-1;
		tmploop=MIN(m1+m2+1,n-k);
		b[i]=0.;
		for (j=MAX(1,1-k); j<=tmploop; j++)
			b[i] += a[i][j]*x[j+k];
	}
}

void bandLU(float **a, int n, int m1, int m2, float **a1, int *index)
/*< band matrix A=LU >*/
{
	int i, j, k, l, mm;
	float dum;

	mm=m1+m2+1;
	l=m1;
	for (i=1; i<=m1; i++){
		for (j=m1+2-i; j<=mm; j++) a[i][j-l]=a[i][j];
		l--;
		for (j=mm-l; j<=mm; j++) a[i][j]=0.;
	}

	l=m1;
	for (k=1; k<=n; k++){
		dum=a[k][1];
		i=k;
		if(l<n) l++;
		for(j=k+1; j<=l; j++){
			if (fabsf(a[j][1]) >fabsf(dum)){
				dum=a[j][1];
				i=j;
			}
		}
		index[k]=i;
		if (dum == 0.) a[k][1]=TINY;
		if(i!=k){
			for(j=1; j<=mm; j++) SWAP(a[k][j],a[i][j])
		}
		for(i=k+1; i<=l; i++){
			dum=a[i][1]/a[k][1];
			a1[k][i-k]=dum;
			for (j=2; j<=mm; j++) a[i][j-1]=a[i][j]-dum*a[k][j];
			a[i][mm]=0.;
		}
	}
}

void bandLUxb(float **a, int n, int m1, int m2, float **a1, int *index, float *b)
/*< band matrix LUx=b >*/
{
	int i,k,l,mm;
	float dum;

	mm=m1+m2+1;
	l=m1;

	for (k=1; k<=n; k++){ /* Forward substitution */
		i=index[k];
		if (i!=k) SWAP(b[k],b[i])
		if (l<n) l++;
		for (i=k+1; i<=l; i++) b[i] -= a1[k][i-k]*b[k];
	}

	l=1;
	for (i=n; i>=1; i--){ /* Backward substitution */
		dum=b[i];
		for (k=2; k<=l; k++) dum -= a[i][k]*b[k+i-1];
		b[i]=dum/a[i][1];
		if(l<mm) l++;
	}
}

void bandSmooth(int n, int n1, float **a)
/*< generate double triangle smoothing operator in band maxtrix form >*/
{
	int i,j,k,l, m1, m2;
	float *p, *dp, sum=0.;

	m1=2*n+1;
	m2=4*n+1;
	p=sf_floatalloc(m1+1);
	dp=sf_floatalloc(m2+1);
	memset(dp, 0., (m2+1)*sizeof(float));

	/* construct filter */
	for(i=1; i<=n+1; i++){
		p[i]=i;
		p[m1+1-i]=i;
	}
	for(i=1; i<=m1; i++){
		for (j=1; j<=m1; j++){
			k=i+j-1;
			dp[k] += p[i]*p[j];
		}
	}

	/* normalize dp */
	for(i=1; i<=m2; i++)
		sum += dp[i];
	for(i=1; i<=m2; i++)
		dp[i] /= sum;

	/* construct band matrix */
	memset(a[0], 0., (n1+1)*(m2+1)*sizeof(float));
	for(i=1; i<=n1; i++){
		for(j=1; j<=m2; j++){
			a[i][j]=dp[j];
		}
	}
	for (i=1; i<=m1-1; i++){
		for(j=1; j<=m1-i; j++){
			k=m2+2-2*i-j;
			a[i][k]=a[i][k]+a[i][j];
			a[i][j]=0.;
		}
	}
	for(i=1; i<=m1-1; i++){
		for(j=1; j<=m2; j++){
			k=n1+1-i;
			l=m2+1-j;
			a[k][l]=a[i][j];
		}
	}

	free(p); free(dp);
}

void bandcp(float *a, int n, float *b)
/*< copy a->b >*/
{
	int i;
	for (i=0; i<n; i++)
		b[i]=a[i];
}

void sf_divn_combine2 (int n, const float* one, const float* two, float *prod)
/*< compute product of two divisions >*/
{
    int i;
    float p;

    for (i=0; i < n; i++) {
	p = sqrtf(fabsf(one[i]*two[i]));
	if ((one[i] > 0. && two[i] < 0. && -two[i] >= one[i]) ||
	    (one[i] < 0. && two[i] > 0. && two[i] >= -one[i])) 
	    p = -p;
	p += 1.;
	p *= p;
	p *= p/16.;
	prod[i] = p;	
    }
}
 
void sf_divn_combine3 (int n, const float* one, const float* two, float *prod)
/*< compute product of two divisions >*/
{
    int i;
    float p;

    for (i=0; i < n; i++) {
	p = one[i]*two[i];
	if(p<0.) p=0.;
	prod[i] = p;	
    }
}

int main(int argc, char* argv[])
{
    bool verb, shift, adjsrc;
	sf_map4 mo;
    int i1, iw, ir, im, k;
    int n1, nw, nr, nd, n12, rect1, m1, mm, tshift, **index;
    float **d, **dd, **p, **pp, **f1, **f2, **f, **g1, **g2,   *str, **smooth, ***a, ***a1, **b, **tmp, meand, meanp, d1, o1, dw, w0, w, dr, r0;
    sf_file dat, flt, mat, adj1, adj2;

    sf_init(argc,argv);

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (!sf_getbool("shift",&shift)) shift=false;
    /* use shift instead of stretch */
    if (!sf_getbool("adjsrc",&adjsrc)) adjsrc=false;
    /* use shift instead of stretch */

    dat = sf_input("obs");
	mat = sf_input("in");
	flt = sf_output("out");
	if(adjsrc){
		adj1 = sf_output("adj1");
		adj2 = sf_output("adj2");
	}

	if(SF_FLOAT != sf_gettype(dat)) sf_error("Need float input");
	if(!sf_histint(dat,"n1",&n1)) sf_error ("No n1= in input");
	if(!sf_histfloat(dat,"d1",&d1)) sf_error ("No d1= in input");
	if(!sf_histfloat(dat,"o1",&o1)) o1 = 0.;

	if(!sf_histint(dat,"n2", &nr)) nr=1; 
	if(!sf_histfloat(dat,"d2", &dr)) dr=1; 
	if(!sf_histfloat(dat,"o2", &r0)) r0=0; 
 
    if (!sf_getint("rect1",&rect1)) rect1=50;
    /* smoothing along first axis */
	if (!sf_getint("nw", &nw)) sf_error("Need nw=");
	/* number of omega values */
	if (!sf_getfloat("dw", &dw)) sf_error("Need dw=");
	/* omega sampling */
	if (!sf_getfloat("w0", &w0)) sf_error("Need w0=");
	/* omega origin */

	nd=n1*nr;
	n12=n1*nw;
	m1=2*rect1;
	mm=4*rect1+1;

	sf_putint(flt, "n2", nw);
	sf_putfloat(flt, "d2", dw);
	sf_putfloat(flt, "o2", w0);
	sf_putstring(flt, "label2", "Gamma");
	sf_putint(flt, "n3", nr);
	sf_putfloat(flt, "d3", dr);
	sf_putfloat(flt, "o3", r0);
	sf_putstring(flt, "label3", "Distance");
	sf_putstring(flt, "unit3", "km");

    d = sf_floatalloc2(n1, nr);
    dd = sf_floatalloc2(n1, nw);
    p = sf_floatalloc2(n1, nr);
    pp = sf_floatalloc2(n1, nw);
    f1 = sf_floatalloc2(n1, nw);
    f2 = sf_floatalloc2(n1, nw);
    f = sf_floatalloc2(n1, nw);
    g1 = sf_floatalloc2(n1, nw);
    g2 = sf_floatalloc2(n1, nw);

    smooth = sf_floatalloc2(mm+1, n1+1);
    a = sf_floatalloc3(mm+1, n1+1, nw);
    a1 = sf_floatalloc3(m1+1, n1+1, nw);
    b = sf_floatalloc2(n1+1, nw);
    tmp = sf_floatalloc2(n1+1, nw);
    index = sf_intalloc2(n1+1, nw);
	bandSmooth(rect1, n1, smooth);

    sf_floatread(d[0],nd,dat);
	sf_floatread(p[0],nd,mat);

	for (ir=0; ir<nr; ir++){
		if(verb) sf_warning("ir/nr=%d/%d", ir+1, nr);

		/* construct reference data and matched data */
		memset(dd[0], 0., n12*sizeof(float));
		if(shift){
			for (iw=0; iw<=nw/2; iw++){
				tshift=(nw/2-iw)*dw/d1;
				for(i1=tshift; i1<n1; i1++)
					dd[iw][i1-tshift]=d[ir][i1];
			}
			for (iw=nw/2+1; iw<nw; iw++){
				tshift=(iw-nw/2)*dw/d1;
				for(i1=0; i1<n1-tshift; i1++)
					dd[iw][i1+tshift]=d[ir][i1];
			}
		}else{
			str = sf_floatalloc(n1);
			mo = sf_stretch4_init(n1, o1, d1, n1, 0.01);
			for (iw=0; iw<nw; iw++){
				for (i1=0; i1<n1; i1++){
					str[i1]=(i1*d1+o1)*(iw*dw+w0);
				}
				sf_stretch4_define(mo, str);
				sf_stretch4_apply(false, mo, d[ir], dd[iw]);
			}
			sf_stretch4_close(mo);
		}

		for (iw=0; iw<nw; iw++){
			for (i1=0; i1<n1; i1++){
				pp[iw][i1]=p[ir][i1];
			}
		}

		/* normalization */
		meand=0.;
		meanp=0.;
		for(iw=0; iw<nw; iw++){
			for (i1=0; i1<n1; i1++){
				meand += dd[iw][i1]*dd[iw][i1];
				meanp += pp[iw][i1]*pp[iw][i1];
			}
		}
		meand = sqrtf(meand/n12);
		meanp = sqrtf(meanp/n12);

#ifdef _OPENMP
#pragma omp parallel for \
		private(iw,i1,im,k) \
		shared(meand,meanp,n1) 
#endif
		for(iw=0; iw<nw; iw++){

			/* 1 */
			/* calculate D(rt) w1 = p(t) */

			/* right-hand side */
			for(i1=0; i1<n1; i1++)
				tmp[iw][i1+1]=dd[iw][i1]*pp[iw][i1]/meand/meand;
			bandcp(smooth[0], (n1+1)*(mm+1), a[iw][0]);
			bandAx(a[iw], n1, m1, m1, tmp[iw], b[iw]);

			/* full matrix */
			for(i1=0; i1<n1; i1++)
				tmp[iw][i1+1]=dd[iw][i1]*dd[iw][i1]/meand/meand-1.;
			for(i1=1; i1<=n1; i1++){
				for(im=1; im<=mm; im++){
					k=m1+1+i1-im;
					if(k>=1 && k<=n1)
						a[iw][k][im] *= tmp[iw][i1];
				}
			}
			for(i1=1; i1<=n1; i1++){
				a[iw][i1][m1+1] = a[iw][i1][m1+1]+1.;
			}

			/* LU decomposition */
			bandLU(a[iw], n1, m1, m1, a1[iw], index[iw]);
			bandLUxb(a[iw], n1, m1, m1, a1[iw], index[iw], b[iw]);
			for(i1=0; i1<n1; i1++)
				f1[iw][i1]=b[iw][i1+1];

			/* 2 */
			/* calculate P(t) w2 = d(rt) */

			/* right-hand side */
			for(i1=0; i1<n1; i1++)
				tmp[iw][i1+1]=dd[iw][i1]*pp[iw][i1]/meanp/meanp;
			bandcp(smooth[0], (n1+1)*(mm+1), a[iw][0]);
			bandAx(a[iw], n1, m1, m1, tmp[iw], b[iw]);

			/* full matrix */
			for(i1=0; i1<n1; i1++)
				tmp[iw][i1+1]=pp[iw][i1]*pp[iw][i1]/meanp/meanp-1.;
			for(i1=1; i1<=n1; i1++){
				for(im=1; im<=mm; im++){
					k=m1+1+i1-im;
					if(k>=1 && k<=n1)
						a[iw][k][im] *= tmp[iw][i1];
				}
			}
			for(i1=1; i1<=n1; i1++){
				a[iw][i1][m1+1] = a[iw][i1][m1+1]+1.;
			}

			/* LU decomposition */
			bandLU(a[iw], n1, m1, m1, a1[iw], index[iw]);
			bandLUxb(a[iw], n1, m1, m1, a1[iw], index[iw], b[iw]);
			for(i1=0; i1<n1; i1++)
				f2[iw][i1]=b[iw][i1+1];

			/* 3 */
			/* combine */
			sf_divn_combine3(n1, f1[iw], f2[iw], f[iw]);

			if(adjsrc){
				w=w0+iw*dw;

				/* 4 */
				/* adjoint source 1 */

				/* right-hand side */
				for(i1=0; i1<n1; i1++)
					b[iw][i1+1]=f[iw][i1]*(w-1.)*(w-1.)*f2[iw][i1];

				/* full matrix */
				for(i1=0; i1<n1; i1++)
					tmp[iw][i1+1]=dd[iw][i1]*dd[iw][i1]/meand/meand-1.;
				bandcp(smooth[0], (n1+1)*(mm+1), a[iw][0]);
				for (i1=1; i1<=n1; i1++)
					for (im=1; im<=mm; im++)
						a[iw][i1][im] *= tmp[iw][i1];
				for(i1=1; i1<=n1; i1++){
					a[iw][i1][m1+1] = a[iw][i1][m1+1]+1.;
				}

				/* LU decomposition */
				bandLU(a[iw], n1, m1, m1, a1[iw], index[iw]);
				bandLUxb(a[iw], n1, m1, m1, a1[iw], index[iw], b[iw]);

				/* generate adjoint source */
				bandAx(smooth, n1, m1, m1, b[iw], tmp[iw]);
				for(i1=0; i1<n1; i1++)
					g1[iw][i1] = tmp[iw][i1+1]*dd[iw][i1]/meand;

				/* 5 */
				/* adjoint source 2 */

				/* right-hand side */
				for(i1=0; i1<n1; i1++)
					b[iw][i1+1]=f[iw][i1]*(w-1.)*(w-1.)*f1[iw][i1];

				/* full matrix */
				for(i1=0; i1<n1; i1++)
					tmp[iw][i1+1]=pp[iw][i1]*pp[iw][i1]/meanp/meanp-1.;
				bandcp(smooth[0], (n1+1)*(mm+1), a[iw][0]);
				for (i1=1; i1<=n1; i1++)
					for (im=1; im<=mm; im++)
						a[iw][i1][im] *= tmp[iw][i1];
				for(i1=1; i1<=n1; i1++){
					a[iw][i1][m1+1] = a[iw][i1][m1+1]+1.;
				}

				/* LU decomposition */
				bandLU(a[iw], n1, m1, m1, a1[iw], index[iw]);
				bandLUxb(a[iw], n1, m1, m1, a1[iw], index[iw], b[iw]);

				/* generate adjoint source */
				bandAx(smooth, n1, m1, m1, b[iw], tmp[iw]);
				for(i1=0; i1<n1; i1++)
					g2[iw][i1] = tmp[iw][i1+1]*pp[iw][i1]/meanp*f2[iw][i1];
			} // end of adjsrc
		} // end of iw
		
		sf_floatwrite(f[0], n12, flt);

		if(adjsrc){
			memset(d[ir], 0., n1*sizeof(float));
			memset(p[ir], 0., n1*sizeof(float));
			for(iw=0; iw<nw; iw++){
				for(i1=0; i1<n1; i1++){
					d[ir][i1] += g1[iw][i1];
					p[ir][i1] += g2[iw][i1];
				}
			}
			sf_floatwrite(d[ir], n1, adj1);
			sf_floatwrite(p[ir], n1, adj2);
		}
	} // end of ir

	exit(0);
}
