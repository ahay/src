/*

Input: T[m][h][t] and its sf-file  m[], Nm h[], Nh 
Output: a(3) three CRS parameters: a(0) m + a(1)m^2 + a(2) h^2

H_11 = 2 Nx Sum(m^2)
H_12 = 2 Nx Sum(m^4)
H_13 = 2 Sum(x^2) Sum(m)

H_22 = 2 Nx Sum(m^4)
H_23 = 2 Sum(x^2) Sum(m^2)

H_33 = 2 Nm Sum(x^2)

da = - grad H^{-1}
*/
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <assert.h>

#define IMIN(X,Z) ( (X < Z) ? X : Z)
#define IMAX(X,Z) ( (X > Z) ? X : Z)

float sum(float * m, int Nm, int k)
{
    int i;
    float s = 0.f;

    if (1 == k) {
	for (i = 0; i < Nm; i++) 
	    s += m[i];
    }
    else if (2 == k) {
	for (i = 0; i < Nm; i++) 
	    s += m[i] * m[i];
    } else { 
	if (3 == k) {
	    for (i = 0; i < Nm; i++) 
		s += m[i] * m[i] * m[i];	
	} else { 
	    assert(4 == k);
	    for (i = 0; i < Nm; i++) 
		s += m[i] * m[i] * m[i] * m[i];
	}
    }


    return s;
}
// H[3][3]
void hessian(float H[3][3], float * m, int Nm, float * h, int Nh)
{
    H[0][0] = 2 * Nh * sum(m, Nm, 2);
    H[0][1] = 2 * Nh * sum(m, Nm, 3);
    H[0][2] = 2 * sum(h, Nh, 2)  * sum(m, Nm, 1);

    H[1][0] = H[0][1];
    H[1][1] = 2 * Nh * sum(m, Nm, 4);
    H[1][2] = 2 * sum(h, Nh, 2)  * sum(m, Nm, 2);

    H[2][0] = H[0][2];
    H[2][1] = H[1][2];
    H[2][2] = 2 * Nm * sum(h, Nh, 4);

    return;
}
// a[3]
float crs_dt(float a[3], float m, float h)
{
    return a[0] * m + a[1] * m * m + a[2] * h * h;
}
// input: T2[h][m] = T^2 - t_0^2
// output: Tcrs[h][m] and err=Sum Sum (T-Tcrs)^2 dm dh
float f_crs(float **dT2, float a[3], float **dT2crs, int Nm, int Nh, float * m, float * h)
{
    double d, sum = 0.f;
    int im, ih;

    for (im = 0; im < Nm; im++) {

	for (ih = 0; ih < Nh; ih++) {
	    
	    dT2crs[im][ih] = crs_dt(a, m[im], h[ih]);

	    d = dT2[im][ih] - dT2crs[im][ih];

	    sum += d * d;
	}
    }
    return sum;
}
// grad[3]
void f_grad_crs(float grad[3], float * m, int Nm, float * h, int Nh, float **dT2, float **dT2crs)
{
    float delta_dT;
    int im, ih;

    grad[0] = grad[1] = grad[2] = 0.f;

    for (im = 0; im < Nm; im++) {

	for (ih = 0; ih < Nh; ih++) {

	    delta_dT = dT2[im][ih] - dT2crs[im][ih];

	    grad[0] += delta_dT * m[im];

	    grad[1] += delta_dT * m[im] * m[im];

	    grad[2] += delta_dT * h[ih] * h[ih];
	}
    }
    grad[0] *= -2.f;
    grad[1] *= -2.f;
    grad[2] *= -2.f;
}
void inverse_mat(float h[3][3], float i[3][3])
{
    int j, k, l;
    float d = 

	-h[0][2]*h[1][1]*h[2][0] + h[0][1]*h[1][2]*h[2][0] + h[0][2]*h[1][0]*h[2][1]

	-h[0][0]*h[1][2]*h[2][1] - h[0][1]*h[1][0]*h[2][2] + h[0][0]*h[1][1]*h[2][2];

    i[0][0] = (- h[1][2] * h[2][1] + h[1][1] * h[2][2]) / d;
    i[0][1] = (  h[0][2] * h[2][1] - h[0][1] * h[2][2]) / d;
    i[0][2] = (- h[0][2] * h[1][1] + h[0][1] * h[1][2]) / d;

    i[1][0] = i[0][1];
    i[1][1] = (- h[0][2] * h[2][0] + h[0][0] * h[2][2]) / d;
    i[1][2] = (  h[0][2] * h[1][0] - h[0][0] * h[1][2]) / d;

    i[2][0] = i[0][2];
    i[2][1] = i[1][2];
    i[2][2] = (- h[0][1] * h[1][0] + h[0][0] * h[1][1]) / d;

    for (l=0;l<3;l++) {
	for (j=0;j < 3; j++) {
	    d = 0;
	    for (k=0;k<3;k++) {
		d += h[l][k]*i[k][j];
	    }
	    if (l == j) 
		assert(fabs(d-1.f) < 1e-6f);
	    else
		assert(fabs(d) < 1e-6f);
	}	
    }

}
float fit_crs_params(float * m, int Nm, float * h, int Nh, float **dT2, float **dT2crs, float a0[3], float a[3])
{

    (void) f_crs(dT2, a0, dT2crs, Nm, Nh, m, h);

    float Hess[3][3], invH[3][3];
    float grad[3];
    float da[3];
    int i, k;   
 
    hessian(Hess, m, Nm, h, Nh);
 
    f_grad_crs(grad, m, Nm, h, Nh, dT2, dT2crs);

    //Hinv = H^-1
    inverse_mat(Hess, invH);
    // da = -H^-1 * grad
    for(i=0; i < 3; i++) {
	
	da[i] = 0.f;

	for(k=0; k < 3; k++) {
	    da[i] -= invH[i][k] * grad[k];	
	}
    }
    a[0] = a0[0] + da[0];
    a[1] = a0[1] + da[1];
    a[2] = a0[2] + da[2];

    return sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
}
void pick_first_picks(float *** Ttraces, 
		      float ** t, 
		      float ot, float dt, int Nm, int Nh, int Nt)
 // t - the times corresponding to 1st local maximum
{
    int it_max, im, ih, it;
    float t_max;
    for (im=0; im < Nm; im++) {
	for (ih = 0; ih < Nh; ih++) {

	    t_max = -1000.f;
	    it_max = -1;

	    for (it = 0; it < Nt; it++) {
		if (Ttraces[im][ih][it] > t_max) {
		    t_max = Ttraces[im][ih][it];
		    it_max = it;
		}
		/*if (t_max > 0.f && Ttraces[im][ih][it] + 1e-6f < t_max) {
		    break;
		    }*/
	    }
	    assert(it_max > 0);
	    t[im][ih] = ot + dt * it_max;
	}
    }
}

int main(int argc, char* argv[])
{
    int  Nm, Nh, im0, Nt, im, ih;
    float dm, om, dh, oh, dt, ot, x0, grad_norm;
    
    float **t;                    /* surface to fit */
    float ***Ttraces, **dT2, **dT2crs, **tcrs;
    float t0, t02;
    float * m_mids, *h_halfoffset;
  
    sf_file in,out_tcrs = NULL,out_tcrs_params = NULL;
    char * out_tcrs_file = 0;
    float crs_a0[3], crs_a[3];
 
    sf_init(argc,argv);

    in = sf_input("in");

 
    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    if (!sf_histint(in,"n3",&Nm)) sf_error("No n3="); //SLOWEST direction, KM
    if (!sf_histfloat(in,"d3",&dm)) sf_error("No d3=");
    if (!sf_histfloat(in,"o3",&om)) sf_error("No o3=");

    if (!sf_histint(in,"n2",&Nh)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2=");

    if (!sf_histint(in,"n1",&Nt)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1=");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1=");
    

    /* memory allocations */
    Ttraces = sf_floatalloc3(Nt,Nh, Nm);

    /* Solution vector */
    t = sf_floatalloc2(Nh,Nm);
    dT2 = sf_floatalloc2(Nh,Nm);
    dT2crs = sf_floatalloc2(Nh,Nm);
    tcrs = sf_floatalloc2(Nh,Nm);

    m_mids = sf_floatalloc(Nm);
    h_halfoffset = sf_floatalloc(Nh);

   
    /* read input escape variable - initial guess */
    sf_floatread(Ttraces[0][0],Nm*Nh*Nt,in);
    
    pick_first_picks(Ttraces, t, ot, dt, Nm,  Nh,  Nt);
    // t - the times corresponding to 1st local maximum
    
    out_tcrs_file = sf_getstring ("tcrs");
    out_tcrs_params = sf_output ("out");
    if (out_tcrs_file) {
	
        out_tcrs = sf_output (out_tcrs_file);

        sf_putint (out_tcrs, "n3", 1);
        sf_putint (out_tcrs, "n2", Nm);
        sf_putint (out_tcrs, "n1", Nh);

        sf_putfloat (out_tcrs, "d3", 0);
        sf_putfloat (out_tcrs, "d2", dm);
        sf_putfloat (out_tcrs, "d1", dh);

        sf_putfloat (out_tcrs, "o3", 0); 
        sf_putfloat (out_tcrs, "o2", om);
        sf_putfloat (out_tcrs, "o1", oh);
    }
    

    /* output */
    //sf_floatwrite(t[0],Nm*Nh,out);

    /* */
    /*m_all_ind = 1 + floor ( 0.5 + (x0 - m_all(1)) / (m_all(2)-m_all(1)) );
      x0 = m_all(m_all_ind);
      t0 = Tmx_all (m_all_ind, 1);*/
    //x0 = m_mids[(int)(Nm/2 + 0.5)];

    im0 = (int)(Nm/2.f + 0.5f) - 1;
    x0 = om + dm * im0;
    t0 = t[im0][0];
    t02 = t0 * t0;
    for (im = 0; im < Nm; im ++) {
	m_mids[im] = /*1000.f* */( (om + im * dm) - x0);
    }

    for (ih = 0; ih < Nh; ih ++) {
	h_halfoffset[ih] = /*1000.f* */(oh + ih*dh);
    }

    for (im = 0; im < Nm; im ++) 
	for (ih = 0; ih < Nh; ih ++) 		
	    dT2[im][ih] = t[im][ih]*t[im][ih] - t02;

    crs_a0[0] = crs_a0[1] = crs_a0[2] = 0.f;

    grad_norm = fit_crs_params(m_mids, Nm, h_halfoffset, Nh, dT2, dT2crs, crs_a0, crs_a);
    sf_warning("fitcrs: INITIAL norm grad = %g crs params for {m, m m, h h} are: [%g %g %g]", grad_norm, crs_a0[0], crs_a0[1], crs_a0[2]);
    grad_norm = fit_crs_params(m_mids, Nm, h_halfoffset, Nh, dT2, dT2crs, crs_a, crs_a0);

    grad_norm = fit_crs_params(m_mids, Nm, h_halfoffset, Nh, dT2, dT2crs, crs_a0, crs_a);
    grad_norm = fit_crs_params(m_mids, Nm, h_halfoffset, Nh, dT2, dT2crs, crs_a, crs_a0);

    sf_warning("fitcrs : HESSIAN MINIMIZATION norm grad = %g crs params for {m, m m, h h} are: [%g %g %g]", grad_norm, crs_a0[0], crs_a0[1], crs_a0[2]);


	sf_floatwrite(crs_a0, 3, out_tcrs_params);
	sf_putint (out_tcrs_params, "n3", 1);
	sf_putint (out_tcrs_params, "n2", 1);
	sf_putint (out_tcrs_params, "n1", 3);

    if (out_tcrs_file) {
	// BUGBUG
	 //sf_intwrite(t_colors[0][0],nz*nx*na,out_colors);
	for (im = 0; im < Nm; im ++) 
	    for (ih = 0; ih < Nh; ih ++) 		
		tcrs[im][ih] = sqrt(dT2crs[im][ih] + t02);

       sf_floatwrite(tcrs[0],Nm*Nh,out_tcrs);

    }

    /* TBD delete !!!
       for (int ia = 0; ia < na; ia++) {
	for (int ix = 0; ix < nx; ix ++) {
	    free(&(t[ia][ix][0]));
	    free(&(t_colors[ia][ix][0]));
	}
	free(t[ia]);
	free(t_colors[ia]);
    }
    */
/*    free(t[0][0]);
    free(t_colors[0][0]);
	
    free(t);
    free(t_colors);

	free(dist_z[0][0]);
	free(dist_z);

	free(dist_x[0][0]);
	free(dist_x);
*/
    exit(0);
}


