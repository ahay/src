/* stiffness tensor for 3D TTI models*/
/*
  Copyright (C) 2007 Colorado School of Mines
  
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
#define COSNU4 pow(cos(n),4)
#define COSNU3 pow(cos(n),3)
#define COSNU2 pow(cos(n),2)

#define SINNU4 pow(sin(n),4)
#define SINNU3 pow(sin(n),3)
#define SINNU2 pow(sin(n),2)

#define COS4NU cos(4*n)
#define COS2NU cos(2*n)
#define COSNU  cos(n)
#define SIN4NU sin(4*n)
#define SIN2NU sin(2*n)
#define SINNU  sin(n)

/*The program output stiffness tensor for TI media in 2D and 3D
  in this order:
  for 2D: C11,C13,C15,
              C33,C35,
	          C55
  for 3d: C11,C12,C13,C14,C15,C16,
              C22,C23,C24,C25,C26,
	          C33,C34,C35,C36,
		      C44,C45,C46,
		          C55,C56,
			      c66

*/
int main(int argc, char* argv[])
{
    bool verb;
    int  dim;

    /* I/O files */
    sf_file Fvp=NULL; 
    sf_file Fvs=NULL; 
    sf_file Fro=NULL; 
    sf_file Fepsilon=NULL; 
    sf_file Fdelta=NULL; 
    sf_file Fgamma=NULL; 
    sf_file Fnu=NULL; 
    sf_file Falpha=NULL; 

    sf_file Fcc=NULL; 
    /* I/O arrays */
    sf_axis az,ax,ay,ac;
    int nc,n1,n2,n3;
    int ix,iy,iz;
    float c33,c55,c11,c13,c66,c12;
    float m11,m12,m13,m15,m22,m23,m25,m33,m35,m44,m46,m55,m66;
    float n;

    float ***vp=NULL;        
    float ***vs=NULL;        
    float ***ro=NULL;        
    float ***epsilon=NULL;        
    float ***delta=NULL;       
    float ***gamma=NULL;     
    float ***nu=NULL;        
    float ***alpha=NULL;   
	
    float ****cc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint("dim",&dim)) dim=2; /* verbosity flag */
    Fvp      = sf_input ("vp" ); 
    Fvs      = sf_input ("vs" ); 
    Fro      = sf_input ("ro" ); 
    Fepsilon = sf_input ("epsilon" ); 
    Fdelta   = sf_input ("delta" ); 
    Fgamma   = sf_input ("gamma" ); 
    Fnu      = sf_input ("nu" ); 
    Falpha   = sf_input ("alpha" ); 
    
    Fcc   = sf_output ("out" ); 
    /* axes */
    az = sf_iaxa(Fvp,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z coordinates*/
    ax = sf_iaxa(Fvp,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x coordinates*/
    ay = sf_iaxa(Fvp,3); sf_setlabel(ay,"y"); /* x coordinates*/
    n1 = sf_n(az);    n2 = sf_n(ax); n3=sf_n(ay);
	
    /* cube axes */

 
    vp=     sf_floatalloc3(n1,n2,n3);  sf_floatread(vp[0][0],     n1*n2*n3, Fvp );
    vs=     sf_floatalloc3(n1,n2,n3);  sf_floatread(vs[0][0],     n1*n2*n3, Fvs );
    ro=     sf_floatalloc3(n1,n2,n3);  sf_floatread(ro[0][0],     n1*n2*n3, Fro );
    epsilon=sf_floatalloc3(n1,n2,n3);  sf_floatread(epsilon[0][0],n1*n2*n3, Fepsilon );
    delta=  sf_floatalloc3(n1,n2,n3);  sf_floatread(delta[0][0],  n1*n2*n3, Fdelta );
    gamma=  sf_floatalloc3(n1,n2,n3);  sf_floatread(gamma[0][0],  n1*n2*n3, Fgamma );
    nu=     sf_floatalloc3(n1,n2,n3);  sf_floatread(nu[0][0],     n1*n2*n3, Fnu );
    alpha=  sf_floatalloc3(n1,n2,n3);  sf_floatread(alpha[0][0],  n1*n2*n3, Falpha );
  
    
    /* setup output data header */

    if(dim==2){
	nc=6; cc=sf_floatalloc4(n1,n2,nc,1);
	ac=sf_maxa(nc,0,1);if(verb) sf_raxa(ac);
	sf_oaxa(Fcc,az,1);sf_oaxa(Fcc,ax,2);sf_oaxa(Fcc,ac,3);
	
	for (ix=0;ix<n2;ix++){
	    for (iz=0;iz<n1;iz++){
		/*2d vti coefficients	    */
		c33=vp[0][ix][iz]*vp[0][ix][iz]*ro[0][ix][iz];
		c55=vs[0][ix][iz]*vs[0][ix][iz]*ro[0][ix][iz];
		c11=2*epsilon[0][ix][iz]*c33+c33;
		c13=sqrt(2*c33*(c33-c55)*delta[0][ix][iz]+(c33-c55)*(c33-c55))-c55;
		/*2d tti coefficients	    */
		n=nu[0][ix][iz]*SF_PI/180.0;
		cc[0][0][ix][iz]=c11*COSNU4+2*(c13+2*c55)*COSNU2*SINNU2+c33*SINNU4;/*c11 of tti*/
		cc[0][1][ix][iz]=(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*COS4NU)/8;/*c13 of tti*/
		cc[0][2][ix][iz]=(c11-c33+(c11-2*c13+c33-4*c55)*COS2NU)*SIN2NU/4;/*c15 of tti*/
		cc[0][3][ix][iz]=c33*COSNU4+2*(c13+2*c55)*COSNU2*SINNU2+c11*SINNU4;/*c33 of tti*/
		cc[0][4][ix][iz]=(c11-c33-(c11-2*c13+c33-4*c55)*COS2NU)*SIN2NU/4;/*c35 of tti*/
		cc[0][5][ix][iz]=(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*COS4NU)/8;/*c55 of tti*/
		
	    }
	}
	sf_floatwrite(cc[0][0][0],n1*n2*nc*1,Fcc);
    }
    else { if(verb) sf_raxa(ay);
	nc=21;cc=sf_floatalloc4(n1,n2,n3,nc);
	ac=sf_maxa(nc,0,1);	if(verb) sf_raxa(ac);
	sf_oaxa(Fcc,az,1);sf_oaxa(Fcc,ax,2);sf_oaxa(Fcc,ay,3);sf_oaxa(Fcc,ac,4);

#ifdef _OPENMP
#pragma omp parallel for						\
    private(ix,iy,iz,c11,c12,c13,c33,c55,c66, \
	    m11,m12,m13,m15,m22,m23,m25,m33,m35,m44,m46,m55,m66,n)	\
    shared(cc)
#endif

	for (iy=0;iy<n3;iy++){	    
	    for (ix=0;ix<n2;ix++){
		for (iz=0;iz<n1;iz++){
		    /*2d vti coefficients	    */
		    c33=vp[iy][ix][iz]*vp[iy][ix][iz]*ro[iy][ix][iz];
		    c55=vs[iy][ix][iz]*vs[iy][ix][iz]*ro[iy][ix][iz];
		    c11=2*epsilon[iy][ix][iz]*c33+c33;
		    c13=sqrt(2*c33*(c33-c55)*delta[iy][ix][iz]+(c33-c55)*(c33-c55))-c55;
		    c66=2*gamma[iy][ix][iz]*c55+c55;
		    c12=c11-2*c66;
		    /*2d tti coefficients	    */
		    n=nu[iy][ix][iz]*SF_PI/180.0;
		    m11=c11*COSNU4+2*(c13+2*c55)*COSNU2*SINNU2+c33*SINNU4;/*c11 of 2d tti*/
		    m12=c12*COSNU2+c13*SINNU2;/*c12 of tti*/
		    m13=(c11+6*c13+c33-4*c55-(c11-2*c13+c33-4*c55)*COS4NU)/8;/*c13 of 2d tti*/
		    m15=(c11-c33+(c11-2*c13+c33-4*c55)*COS2NU)*SIN2NU/4;/*c15 of 2d tti*/
		    m22=c11;/*c22 of 2d tti*/
		    m23=c13*COSNU2+c12*SINNU2;
		    m25=(c12-c13)*COSNU*SINNU;
		    m33=c33*COSNU4+2*(c13+2*c55)*COSNU2*SINNU2+c11*SINNU4;/*c33 of 2d tti*/
		    m35=-(-c11+c33+(c11-2*c13+c33-4*c55)*COS2NU)*SIN2NU/4;/*c35 of 2d tti*/
		    m44=c55*COSNU2+c66*SINNU2;/*c44 of 2d tti*/
		    m46=(c66-c55)*COSNU*SINNU;/*c46 of tti*/
		    m55=(c11-2*c13+c33+4*c55-(c11-2*c13+c33-4*c55)*COS4NU)/8;/*c55 of tti*/
		    m66=c66*COSNU2+c55*SINNU2;/*c66 of tti*/

		    n=alpha[iy][ix][iz]*SF_PI/180.0;
		    cc[0][iy][ix][iz]=m11*COSNU4+2*(m12+2*m66)*COSNU2*SINNU2+m22*SINNU4;/*c11 of 3d tti*/
		    cc[1][iy][ix][iz]=(m11+6*m12+m22-4*m66-(m11-2*m12+m22-4*m66)*COSNU4)/8;/*c12 of 3d tti*/
		    cc[2][iy][ix][iz]=m13*COSNU2+m23*SINNU2;/*c13 of 3d tti*/
		    cc[3][iy][ix][iz]=(m15-2*m46)*COSNU2*SINNU+m25*SINNU3;/*c14 of 3d tti*/
		    cc[4][iy][ix][iz]=m15*COSNU3+(m25+2*m46)*COSNU*SINNU2;/*c15 of 3d tti*/
		    cc[5][iy][ix][iz]=(m11-m22-(m11-2*m12+m22-4*m66)*COS2NU)*SIN2NU/4;/*c16 of 3d tti*/
		    cc[6][iy][ix][iz]=m22*COSNU4+2*(m12+2*m66)*COSNU2*SINNU2+m11*SINNU4;/*c22 of 3d tti*/
		    cc[7][iy][ix][iz]=m23*COSNU2+m13*SINNU2;/*c23 of 3d tti*/
		    cc[8][iy][ix][iz]=(m25+2*m46)*COSNU2*SINNU+m15*SINNU3;/*c24 of 3d tti*/
		    cc[9][iy][ix][iz]=m25*COSNU3+(m15-2*m46)*COSNU*SINNU2;/*c25 of 3d tti*/
		    cc[10][iy][ix][iz]=(m11-m22-(m11-2*m12+m22-4*m66)*COS2NU)*SIN2NU/4;/*c26 of 3d tti*/
		    cc[11][iy][ix][iz]=m33;/*c33 of 3d tti*/
		    cc[12][iy][ix][iz]=m35*SINNU;/*c34 of 3d tti*/
		    cc[13][iy][ix][iz]=m35*COSNU;/*c35 of 3d tti*/
		    cc[14][iy][ix][iz]=(m13-m23)*COSNU*SINNU;/*c36 of 3d tti*/
		    cc[15][iy][ix][iz]=m44*COSNU2+m55*SINNU2;/*c44 of 3d tti*/
		    cc[16][iy][ix][iz]=(-m44+m55)*COSNU*SINNU;/*c45 of 3d tti*/
		    cc[17][iy][ix][iz]=m46*COSNU3+(m15-m25-m46)*COSNU*SINNU2;/*c46 of 3d tti*/
		    cc[18][iy][ix][iz]=m55*COSNU2+m44*SINNU2;/*c55 of 3d tti*/
		    cc[19][iy][ix][iz]=m46*SINNU3+(m15-m25-m46)*COSNU2*SINNU;/*c56 of 3d tti*/
		    cc[20][iy][ix][iz]=(m11-2*m12+m22+4*m66-(m11-2*m12+m22-4*m66)*COS4NU)/8;/*c66 of 3d tti*/

		}
	    }
	}
	sf_floatwrite(cc[0][0][0],n1*n2*n3*nc,Fcc);
    }
  
    
    exit (0);
}
