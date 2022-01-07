/* Calculate CRE trajectory on CMP x Offset plane given zero-offset CRS parameters (RN, RNIP, BETA)

Programmer: Rodolfo A. C. Neves (Dirack) 31/08/2019

Email:  rodolfo_profissional@hotmail.com 

License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>. 

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>

int main(int argc, char* argv[])
{

	float** m; // CMP
	float h; // Half offset
	float m0; // Central CMP
	double alpha; // Assymetry parameter
	float om0; // m0's axis origin
	float dm0; // m0's sampling
	int nm0; // Number of m0s
	float ot0; // t0's axis origin
	float dt0; // t0's sampling
	int nt0; // Number of t0's
	float** p; // RNIP and BETA parameters temporary vector
	int np1; // Number of parameters in parameters file
	int np2; // Number of (t0, m0) pairs in parameters file
	bool verb; // Verbose parameter
	float dm; // CMP sampling
	float om; // CMP axis origin
	int nm; // Number of CMP samples
	float dh; // Half offset sampling
	float oh; // Half offset axis origin
	int nh; // Number of Half offset samples
	float dt; // Time sampling
	float ot; // Time axis origin
	int nt; // Number of time samples
	int i,l,k; // counters

	/* RSF files I/O */  
	sf_file in, out, par;

	/* RSF files axis */
	sf_axis ax,ay,az;

	sf_init(argc,argv);

	in = sf_input("in"); // Data cube A(m,h,t)
	par = sf_input("param"); // RNIP and BETA parameters
	out = sf_output("out"); // m(h) vector CRE coordinates

	if (!sf_getint("nm0",&nm0)) sf_error("Need nm0");
	/* Number of central CMPs in parameters file */
	if (!sf_getfloat("om0",&om0)) sf_error("Need om0");
	/* First central CMP coordinate in parameters file (Km) */
	if (!sf_getfloat("dm0",&dm0)) sf_error("Need dm0");
	/* Central CMPs sampling in parameters file (Km) */
	if (!sf_getint("nt0",&nt0)) sf_error("Need nt0");
	/* Number of t0s in parameters file */
	if (!sf_getfloat("ot0",&ot0)) sf_error("Need ot0");
	/* First t0 coordinate in parameters file (s) */
	if (!sf_getfloat("dt0",&dt0)) sf_error("Need dt0"); 
	/* t0s sampling in parameters file (s) */

	/* Parameters file */
	if (!sf_histint(par,"n1",&np1)) sf_error("No n1= in parameters file");
	if (!sf_histint(par,"n2",&np2)) sf_error("No n2= in parameters file");

	/* Check dimensions and input */
	if ((nm0*nt0) != np2) 
		sf_error("nm0*nt0 should be equal to n2 in parameters file!");
	if(np1 < 3) 
		sf_error("It should have at least 3 parameters in parameters file!");

	/* seismic data cube A(m,h,t) */
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input file");
	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input file");
	if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input file");
	if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input file");
	if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input file");
	if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input file");
	if (!sf_histint(in,"n3",&nm)) sf_error("No n3= in input file");
	if (!sf_histfloat(in,"d3",&dm)) sf_error("No d3= in input file");
	if (!sf_histfloat(in,"o3",&om)) sf_error("No o3= in input file");

	if(! sf_getbool("verb",&verb)) verb=0;
	/* 1: active mode; 0: quiet mode */

	if (verb) {

		sf_warning("Active mode on!!!");
		sf_warning("Command line parameters: "); 
		sf_warning("nm0=%d om0=%f dm0=%f",nm0,om0,dm0);
		sf_warning("nt0=%d ot0=%f dt0=%f",nt0,ot0,dt0);
		sf_warning("Parameters file: ");
		sf_warning("n1=%d n2=%d",np1,np2);
		sf_warning("Data cube dimensions: ");
		sf_warning("n1=%d d1=%f o1=%f",nt,dt,ot);
		sf_warning("n2=%d d2=%f o2=%f",nh,dh,oh);
		sf_warning("n3=%d d3=%f o3=%f",nm,dm,om);
	}
	
	m = sf_floatalloc2(nh,np2);
	p = sf_floatalloc2(np1,np2);
	sf_floatread(p[0],np1*np2,par);

	for(l=0;l<nm0;l++){

		m0 = l*dm0+om0;

		for(k=0;k<nt0;k++){

			alpha = sin(p[(l*nt0)+k][2])/p[(l*nt0)+k][1];

			if(alpha <= 0.001 && alpha >= -0.001){
				for(i=0;i<nh;i++){
					m[(l*nt0)+k][i] = m0;
				}
			}else{
				for(i=0;i<nh;i++){
					h = (dh*i) + oh;
					m[(l*nt0)+k][i] = m0 + (1/(2*alpha)) * (1 - sqrt(1 + 4 * alpha * alpha * h * h));
				}
			}
		}/* loop over t0s */
	}/*loop over m0s */

	/* axis = sf_maxa(n,o,d)*/
	ax = sf_maxa(nh, oh, dh);
	ay = sf_maxa(np2, 0, 1);
	az = sf_maxa(1, 0, 1);

	/* sf_oaxa(file, axis, axis index) */
	sf_oaxa(out,ax,1);
	sf_oaxa(out,ay,2);
	sf_oaxa(out,az,3);
	sf_floatwrite(m[0],nh*np2,out);

	exit(0);
}
