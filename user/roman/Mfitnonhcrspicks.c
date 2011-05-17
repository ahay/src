/*
Compute fitting of Non-hyperbolic CRS to first-arrivals T[m][h].

Input: T[m][h] arrivals and in its sf-file  m[], Nm h[], Nh 

parameters:

A1
A2
B

FF = F(t0,t02,A1,A2,m-h)*F(t0,t02,A1,A2,m+h);

t2 = 0.5 * (F(t0, t02, A1, A2, m) + (2*B+A1*A1-A2)*h*h + sqrt(FF) );

Output: 
     tcrs=filname - t[m][h] crs surface
     prints three CRS parameters a[3] where: t[m][h] = a(0) m + a(1)m^2 + a(2) h^2

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

float F(float t0, float t02, float A1, float A2, float m)
{
    //%F= t02 + 2*A1*t0*m + A2*m*m;
    return (t0 + A1*m)*(t0+A1*m) + A2*m*m;
}

float f_nonhcrs(float m, float h, float A1, float A2, float B, float t0, float t02)
{
    const float FF = F(t0,t02,A1,A2,m-h)*F(t0,t02,A1,A2,m+h);

    if (FF < 0)    
	assert(0);//FF = 0;%-FF;    
    //%t2 = 0.5 * ( F(t0,t02,A1,A2,m) + C + sqrt(FF) );

    return /*t2*/ 0.5 * (F(t0, t02, A1, A2, m) + (2*B+A1*A1-A2)*h*h + sqrt(FF) );
}

float f_hyper_t2(float H2, float sa2, float ca2, float xc, float V_0_2, float m, float x)
{
    float xs = m - x -xc;
    float xs2 = xs*xs;
    
    float xr = m + x - xc;
    float xr2 = xr*xr;
    
    float f = (H2 + xs2*sa2)*(H2 + xr2*sa2);
    /*if (f < 0.f) {
	f
	  V_0_2
	  f = 0;
    end
    }*/
    float t2 = 2*H2 + xs2 + xr2 - 2*xs*xr*ca2+  2*sqrt(f);

    t2 = t2/V_0_2;
    return t2;
}
int main(int argc, char* argv[])
{
    float H2, sa2,  ca2,  xc, alpha, sa, ca, H, V2;
    int  Nm, Nh, im0, im, ih;
    float dm, om, dh, oh, x0,  A1, A2, B, sum, max_err;
    
    float **t;                    /* surface to fit */
    float **dT2, **dT2crs, **tcrs;
    float t0, t02;
    float * m_mids, *h_halfoffset;
  
    sf_file in,/*out,*/out_tcrs;
    //char * out_tcrs_file = 0;
    //double crs_a0[3], crs_a[3], sum;
 
    sf_init(argc,argv);

    in = sf_input("in");
    out_tcrs = sf_output("out");

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    if (!sf_histint(in,"n3",&Nm)) sf_error("No n3="); //SLOWEST direction, KM
    if (!sf_histfloat(in,"d3",&dm)) sf_error("No d3=");
    if (!sf_histfloat(in,"o3",&om)) sf_error("No o3=");
    assert(1 == Nm);

    if (!sf_histint(in,"n2",&Nm)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dm)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&om)) sf_error("No o2=");

    if (!sf_histint(in,"n1",&Nh)) sf_error("No n1=");
    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1=");
    if (!sf_histfloat(in,"o1",&oh)) sf_error("No o1=");
    

    if (!sf_getfloat("A1",&A1)) sf_error("Need A1 (float) A1=");
    if (!sf_getfloat("A2",&A2)) sf_error("Need A2 (float) A2=");
    if (!sf_getfloat("B",&B)) sf_error("Need B (float) B=");

    /* memory allocations */

    /* Solution vector */
    t = sf_floatalloc2(Nh,Nm);
    dT2 = sf_floatalloc2(Nh,Nm);
    dT2crs = sf_floatalloc2(Nh,Nm);
    tcrs = sf_floatalloc2(Nh,Nm);

    m_mids = sf_floatalloc(Nm);
    h_halfoffset = sf_floatalloc(Nh);

   
    /* read input escape variable - initial guess */
    sf_floatread(t[0],Nm*Nh,in);
    
    // t - the times corresponding to 1st local maximum
        
    //if (out_tcrs_file) 
    {
	//out_tcrs = sf_output (out_tcrs_file);

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

    for (im = 0; im < Nm; im ++) {
	for (ih = 0; ih < Nh; ih ++) {
	    dT2[im][ih] = t[im][ih]*t[im][ih] - t02;
	}
    }


    /* output colorsc */
    //if (out_tcrs_file) {
    // BUGBUG
    //sf_intwrite(t_colors[0][0],nz*nx*na,out_colors);
    sum = 0.f;
    max_err = 0.f;

    /*
alpha = a(1);

sa = sin(alpha); 
ca = cos(alpha);
sa2 = sa*sa;
ca2 = ca*ca;

xc = a(2);

%H2 = a(3);
H = a(3);
H2 = H*H;

V2 = 4 * (H2 + xc*xc*sa2) / t02; 

% check: 
%t02 == f_hyper_t2( H2, sa2, ca2, xc, V2, 0, 0);
for im=1:Nm
    for ix=1:Nx
        dT2(im,ix)=f_hyper_t2( H2, sa2, ca2, xc, V2, m(im), x(ix) ) - t02;
    end
end
    */
    alpha = A1;

    sa = sinf(alpha); 
    ca = cosf(alpha);
    sa2 = sa*sa;
    ca2 = ca*ca;

    xc = A2;

    H = B;
    H2 = H*H;

    V2 = 4. * (H2 + xc*xc*sa2) / t02; 

    //check: 
    //t02 == f_hyper_t2( H2, sa2, ca2, xc, V2, 0, 0);
    //dT2(im,ix)=f_hyper_t2( H2, sa2, ca2, xc, V2, m(im), x(ix) ) - t02;

    for (im = 0; im < Nm; im ++) {
	for (ih = 0; ih < Nh; ih ++) {


	    dT2crs[im][ih] = //f_nonhcrs(m_mids[im], h_halfoffset[ih], A1, A2, B, t0, t02);
		f_hyper_t2(H2, sa2,  ca2,  xc,  V2,  
			   m_mids[im], h_halfoffset[ih]);
			   //m,  x);

	    tcrs[im][ih] = sqrt(dT2crs[im][ih]);
	    
	    sum += (tcrs[im][ih] - t[im][ih]) * (tcrs[im][ih] - t[im][ih]);
	    if (max_err < fabs(tcrs[im][ih] - t[im][ih]))
		max_err = fabs(tcrs[im][ih] - t[im][ih]);		
	}
    }
    sf_warning(" Nonhyperbolic crs error-L2 = %g error-max =%g", sqrt(sum/Nm/Nh), max_err);

    sf_floatwrite(tcrs[0],Nm*Nh,out_tcrs);
    
    //}

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


