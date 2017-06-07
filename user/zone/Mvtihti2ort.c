/* Combining VTI and HTI parameters to orthorhombic according to Schoenberg & Sayer (1995)
NOTE: HTI is defined in VTI grid with respect to the vertical (Ruger(1997)) 
Refer to SEAM 2 notes for detailed description

*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>
#include <rsf.h>

/* definition for LAPACK ROUTINE */
int     N=6;       
int     LDA=6;     
int     LWORK=360;  
int     INFO;
int     IPIV[6];
double  WORK[360];  

int main(int argc, char* argv[])
{
	int n[3],i,j,nm;
	bool rotate;
	float o[3],d[3];
    float C11v, C12v, C13v, C33v, C55v, C66v, C11h, C12h, C13h, C22h, C44h, C55h, Phi, cortf[36];
    double cvti[36],chti[36], ciso[36], cort[36], bond[36], temp[36];
    sf_file c11o,c22o,c33o,c44o,c55o,c66o,c12o,c13o,c23o,c16o,c26o,c36o,c45o,c33,c55,c11v,c66v,c12v,c13v,c11h,c55h,c12h,c13h,phi;

    sf_init (argc,argv);
    c11o = sf_output("out");
    c22o = sf_output("c22o");
    c33o = sf_output("c33o");
    c44o = sf_output("c44o");
    c55o = sf_output("c55o");
    c66o = sf_output("c66o");
    c12o = sf_output("c12o");
    c13o = sf_output("c13o");
    c23o = sf_output("c23o");

    c33 = sf_input("in");
    c55 = sf_input("c55");
    c11v = sf_input("c11v");
    c66v = sf_input("c66v");
    c12v = sf_input("c12v");
    c13v = sf_input("c13v");
    c11h = sf_input("c11h");
    c55h = sf_input("c55h");
    c12h = sf_input("c12h");
    c13h = sf_input("c13h");
    
    if (!sf_getbool("rotate",&rotate)) rotate=false;
    /* Doing azimuthal rotation (y-> mono, n-> ortho)*/
    if (rotate) {
    	phi = sf_input("phi");
		c16o = sf_output("c16o");
		c26o = sf_output("c26o");
		c36o = sf_output("c36o");
		c45o = sf_output("c45o");
    }	

    /* get 3-D grid parameters */
    if (!sf_histint(c33,"n1",n))     sf_error("No n1= in input");
    if (!sf_histint(c33,"n2",n+1))   n[1]=1;
    if (!sf_histint(c33,"n3",n+2))   n[2]=1;
    if (!sf_histfloat(c33,"d1",d))   sf_error("No d1= in input");
    if (!sf_histfloat(c33,"d2",d+1)) d[1]=1.;
    if (!sf_histfloat(c33,"d3",d+2)) d[2]=1.;
    if (!sf_histfloat(c33,"o1",o))   o[0]=0.;
    if (!sf_histfloat(c33,"o2",o+1)) o[1]=0.;
    if (!sf_histfloat(c33,"o3",o+2)) o[2]=0.;

    /* get all cij (default = orthorhombic) ---------------------------------------------------------*/
    nm = n[0]*n[1]*n[2];

	for (i=0;i<nm;i++) {

	// Read VTI
    sf_floatread(&C33v,1,c33); sf_floatread(&C55v,1,c55); sf_floatread(&C11v,1,c11v); 
    sf_floatread(&C12v,1,c12v); sf_floatread(&C66v,1,c66v); sf_floatread(&C13v,1,c13v); 
    
    // Read HTI
    C22h=C33v; C44h=C55v; sf_floatread(&C11h,1,c11h);  sf_floatread(&C12h,1,c12h); 
    sf_floatread(&C55h,1,c55h); sf_floatread(&C13h,1,c13h);  
    if (rotate) sf_floatread(&Phi,1,phi); 


	// Construct 6x6 array
	cvti[0]= C11v; cvti[6]= C12v;  cvti[12]= C13v; cvti[18]= 0.0 ; cvti[24]= 0.0 ; cvti[30]= 0.0 ; 
	cvti[1]= C12v; cvti[7]= C11v;  cvti[13]= C13v; cvti[19]= 0.0 ; cvti[25]= 0.0 ; cvti[31]= 0.0 ;
	cvti[2]= C13v; cvti[8]= C13v;  cvti[14]= C33v; cvti[20]= 0.0 ; cvti[26]= 0.0 ; cvti[32]= 0.0 ; 
	cvti[3]= 0.0 ; cvti[9]= 0.0 ;  cvti[15]= 0.0 ; cvti[21]= C55v; cvti[27]= 0.0 ; cvti[33]= 0.0 ;
	cvti[4]= 0.0 ; cvti[10]= 0.0 ; cvti[16]= 0.0 ; cvti[22]= 0.0 ; cvti[28]= C55v; cvti[34]= 0.0 ; 
	cvti[5]= 0.0 ; cvti[11]= 0.0 ; cvti[17]= 0.0 ; cvti[23]= 0.0 ; cvti[29]= 0.0 ; cvti[35]= C66v;
	
	chti[0]= C11h; chti[6]= C12h;  chti[12]= C13h; chti[18]= 0.0 ; chti[24]= 0.0 ; chti[30]= 0.0 ; 
	chti[1]= C12h; chti[7]= C22h;  chti[13]= C12h; chti[19]= 0.0 ; chti[25]= 0.0 ; chti[31]= 0.0 ;
	chti[2]= C13h; chti[8]= C12h;  chti[14]= C11h; chti[20]= 0.0 ; chti[26]= 0.0 ; chti[32]= 0.0 ; 
	chti[3]= 0.0 ; chti[9]= 0.0 ;  chti[15]= 0.0 ; chti[21]= C44h; chti[27]= 0.0 ; chti[33]= 0.0 ;
	chti[4]= 0.0 ; chti[10]= 0.0 ; chti[16]= 0.0 ; chti[22]= 0.0 ; chti[28]= C55h; chti[34]= 0.0 ; 
	chti[5]= 0.0 ; chti[11]= 0.0 ; chti[17]= 0.0 ; chti[23]= 0.0 ; chti[29]= 0.0 ; chti[35]= C44h; 
	
	ciso[0]= C33v; ciso[6]= C33v-2*C55v;  ciso[12]= C33v-2*C55v; ciso[18]= 0.0 ; ciso[24]= 0.0 ; ciso[30]= 0.0 ; 
	ciso[1]= C33v-2*C55v; ciso[7]= C33v;  ciso[13]= C33v-2*C55v; ciso[19]= 0.0 ; ciso[25]= 0.0 ; ciso[31]= 0.0 ;
	ciso[2]= C33v-2*C55v; ciso[8]= C33v-2*C55v;  ciso[14]= C33v; ciso[20]= 0.0 ; ciso[26]= 0.0 ; ciso[32]= 0.0 ; 
	ciso[3]= 0.0 ; ciso[9]= 0.0 ;  ciso[15]= 0.0 ; ciso[21]= C55v; ciso[27]= 0.0 ; ciso[33]= 0.0 ;
	ciso[4]= 0.0 ; ciso[10]= 0.0 ; ciso[16]= 0.0 ; ciso[22]= 0.0 ; ciso[28]= C55v; ciso[34]= 0.0 ; 
	ciso[5]= 0.0 ; ciso[11]= 0.0 ; ciso[17]= 0.0 ; ciso[23]= 0.0 ; ciso[29]= 0.0 ; ciso[35]= C55v; 

	// LAPACK inversion
	//VTI
	dgetrf_(&N,&N, cvti, &LDA, IPIV, &INFO);	// LU-factorization
	dgetri_( &N, cvti, &LDA, IPIV, WORK, &LWORK, &INFO ); // General matrix inversion

	//HTI
	dgetrf_(&N,&N, chti, &LDA, IPIV, &INFO);	// LU-factorization
	dgetri_( &N, chti, &LDA, IPIV, WORK, &LWORK, &INFO ); // General matrix inversion
	
	//ISO
	dgetrf_(&N,&N, ciso, &LDA, IPIV, &INFO);	// LU-factorization
	dgetri_( &N, ciso, &LDA, IPIV, WORK, &LWORK, &INFO ); // General matrix inversion
	
	// Orthorhombic compliance
	for (j=0;j<36;j++) {
		cort[j] = cvti[j] + chti[j] - ciso[j];
	}
	
	//ORT (compliance -> stiffness)
	dgetrf_(&N,&N, cort, &LDA, IPIV, &INFO);	// LU-factorization
	dgetri_( &N, cort, &LDA, IPIV, WORK, &LWORK, &INFO ); // General matrix inversion


	 if (rotate) {
		// Bond transformation of resultant orthorhombic (Map view CCW positive w/ Phi wrt x-axis or y-axis)
		double c=cos(Phi), s=sin(Phi), cc=cos(2*Phi), ss=sin(2*Phi);
		bond[0]= c*c ; bond[6]= s*s ;  bond[12]= 0.0 ; bond[18]= 0.0 ; bond[24]= 0.0 ; bond[30]=-ss  ; 
		bond[1]= s*s ; bond[7]= c*c ;  bond[13]= 0.0 ; bond[19]= 0.0 ; bond[25]= 0.0 ; bond[31]= ss  ;
		bond[2]= 0.0 ; bond[8]= 0.0 ;  bond[14]= 1.0 ; bond[20]= 0.0 ; bond[26]= 0.0 ; bond[32]= 0.0 ; 
		bond[3]= 0.0 ; bond[9]= 0.0 ;  bond[15]= 0.0 ; bond[21]= c   ; bond[27]= s   ; bond[33]= 0.0 ;
		bond[4]= 0.0 ; bond[10]= 0.0 ; bond[16]= 0.0 ; bond[22]= -s  ; bond[28]= c   ; bond[34]= 0.0 ; 
		bond[5]= ss/2; bond[11]=-ss/2; bond[17]= 0.0 ; bond[23]= 0.0 ; bond[29]= 0.0 ; bond[35]= cc  ; 

		// C.B^T
		char str[2] = "nt"; double zero=0.0, one=1.0;
		dgemm_(str,str+1,&N,&N,&N,&one,cort,&LDA,bond,&LDA,&zero,temp,&LDA);
	
		// B(C.B^T)
		dgemm_(str,str,&N,&N,&N,&one,bond,&LDA,temp,&LDA,&zero,cort,&LDA);
	}
	
	// Cast to float
	for (j=0;j<36;j++) {
		cortf[j] = (float) cort[j];
	}
	
	sf_floatwrite(cortf,1,c11o); sf_floatwrite(cortf+1,1,c12o); sf_floatwrite(cortf+2,1,c13o);
	sf_floatwrite(cortf+7,1,c22o); sf_floatwrite(cortf+8,1,c23o); sf_floatwrite(cortf+14,1,c33o);
	sf_floatwrite(cortf+21,1,c44o); sf_floatwrite(cortf+28,1,c55o); sf_floatwrite(cortf+35,1,c66o); 
	
	if (rotate) { // Monoclinic
		sf_floatwrite(cortf+5,1,c16o); sf_floatwrite(cortf+11,1,c26o); sf_floatwrite(cortf+17,1,c36o);
		sf_floatwrite(cortf+22,1,c45o);
	}	
	
	}
	
    exit (0);
}

/*         $Id$         */
