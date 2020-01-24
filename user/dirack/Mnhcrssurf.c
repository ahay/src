/* Version 1.0 - Build Non-Hyperbolic CRS approximation surface giver RN, RNIP and BETA parameters.

Programer: Rodolfo A. C. Neves (Dirack) 19/09/2019

Email:  rodolfo_profissional@hotmail.com

License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>

int main(int argc, char* argv[])
{

	float m0; // central CMP
	float om; // CMP axis origin
	float dm; // CMP sampling
	int nm; // Number of CMP's
	float oh; // Offset axis origin
	float dh; // Offset sampling
	int nh; // Number of Offsets
	bool verb; // Key to turn On/Off active mode
	float v0; // Near surface velocity
	float t0; // Normal ray time travel
	float *c; // Temporary parameters vector - last iteration
	float **t; // CRS surface t(m,h)
	float RN, RNIP, BETA; // CRS parameters
	float Fd, Fd1, Fd2;
	float c1, a1, a2, b2;
	float m;
	float h;
	int nc, ih, im;

	/* RSF files I/O */  
	sf_file in, out, par;

	/* RSF files axis */
	sf_axis ax,ay,az;

	sf_init(argc,argv); 

	in = sf_input("in");
	par = sf_input("param");
	out = sf_output("out");

	if (!sf_getfloat("m0",&m0)) m0=0;
	/* central CMP of the approximation (Km) */

	if (!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	if (!sf_getfloat("t0",&t0)) t0=1.5;
	/* Normal ray traveltime (s) */

	if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o1",&oh)) sf_error("No o1= in input");
	if (!sf_histint(in,"n2",&nm)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dm)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&om)) sf_error("No o2= in input");

	if(!sf_histint(par,"n1",&nc)) sf_error("No n1= in parameters input");

	if(! sf_getbool("verb",&verb)) verb=0;
	/* 1: active mode; 0: quiet mode */

	if (verb) {

		sf_warning("Active mode on!!!");
		sf_warning("Command line parameters: "); 
		sf_warning("m0=%f v0=%f t0=%f",m0,v0,t0);
		sf_warning("Input file parameters: ");
		sf_warning("n1=%i d1=%f o1=%f",nh,dh,oh);
		sf_warning("n2=%i d2=%f o2=%f",nm,dm,om);
		sf_warning("Param file parameters: ");
		sf_warning("n1=%i",nc);
	}

	c = sf_floatalloc(nc);
	sf_floatread(c,nc,par);

	RN = c[0];
	RNIP = c[1];
	BETA  = c[2];


	t = sf_floatalloc2(nh,nm);
	
	for (im=0; im < nm; im++){
			
		for(ih=0;ih<nh;ih++){

		m = om + (im * dm);
		m = m - m0;
		h = oh + (ih * dh);
	
		a1=(2*sin(BETA))/(v0);		
		a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RN);
		b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		c1=2*b2+a1*a1-a2;
													
		Fd=(t0+a1*m)*(t0+a1*m)+a2*m*m;				
		Fd2=(t0+a1*(m-h))*(t0+a1*(m-h))+a2*(m-h)*(m-h);
		Fd1=(t0+a1*(m+h))*(t0+a1*(m+h))+a2*(m+h)*(m+h);					
		t[im][ih]=sqrt((Fd+c1*h*h+sqrt(Fd2*Fd1))*0.5); 

		}
	}

	/* axis = sf_maxa(n,o,d)*/
	ax = sf_maxa(nh, oh, dh);
	ay = sf_maxa(nm, om, dm);
	az = sf_maxa(1, 0, 1);

	/* sf_oaxa(file, axis, axis index) */
	sf_oaxa(out,ax,1);
	sf_oaxa(out,ay,2);
	sf_oaxa(out,az,3);
	sf_floatwrite(t[0],nh*nm,out);

	exit(0);
}
