/* Build CRE gather given CMP X Offset CRE trajectory coordinates and interpolated data cube

This program searches for the closest trace to the CRE trajectory to build the CRE Gather for each (m, h) pair given in the interpolated data cube. 

Programmer: Rodolfo A. C. Neves (Dirack) 04/09/2019

Email:  rodolfo_profissional@hotmail.com  

License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>

int main(int argc, char* argv[])
{

	float*** t; // Data cube A(m,h,t)
	float*** creGather; // CRE gathers
	float** m; // CMPs
	bool verb; // Verbose parameters
	float dm; // CMP sampling
	float om; // CMP axis origin
	int nm; // Number of CMP samples
	float dh; // Offset sampling
	float oh; // Offset axis origin
	int nh; // Offset number of samples
	float dt; // Time sampling
	float ot; // Time axis origin
	int nt; // Number of time samples
	int i,j,k,l; // loop counter
	int cn1; // Number of traces in CRE vector
	float cd1; // CRE Gather sampling
	float co1; // CRE Gather axis origin
	int cn2; // Number of m0s x t0s pairs
	int trac_m; // CMP sample index
	int aperture; // Number of traces in CRE Gather
	float mMax; // maximum CMP coordinate of the model
	int nm0; // Number of m0s
	int nt0; // Number of t0s

	/* RSF files I/O */  
	sf_file in, out, out_m, cremh;

	/* RSF files axis */
	sf_axis ax,ay,az,am1,am2;

	sf_init(argc,argv);

	in = sf_input("in"); // Data cube A(m,h,t)
	cremh = sf_input("cremh"); /* CRE Gather m(h) coordinates */
	out = sf_output("out"); // CRE Gather
	out_m = sf_output("m"); // CRE Gather CMP coordinate

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

	/* cre trajectories m(h) */
	if(!sf_histint(cremh,"n1",&cn1)) sf_error("No n1= in cremh file");
	if(!sf_histfloat(cremh,"d1",&cd1)) sf_error("No d1= in cremh file");
	if(!sf_histfloat(cremh,"o1",&co1)) sf_error("No o1= in cremh file");
	if(!sf_histint(cremh,"n2",&cn2)) sf_error("No n2= in cremh file");

	if(!sf_getint("nm0",&nm0)) sf_error("Need nm0");
	/* Number of central CMPs in cremh file */
	if(!sf_getint("nt0",&nt0)) sf_error("Need nt0");
	/* Number of t0s in cremh file */
	if(!sf_getint("aperture",&aperture)) aperture=1;
	/* Number of traces to put in a CRE Gather*/

	if(aperture > cn1){
		sf_error("The aperture can't be > n1 in cremh file\naperture=%i n2=%i",aperture,cn1);
	}

	if(! sf_getbool("verb",&verb)) verb=0;
	/* 1: active mode; 0: quiet mode */

	if (verb) {

		sf_warning("Active mode on!!!");
		sf_warning("CRE gather coordinates m(h) (cremh file): ");
		sf_warning("n1=%d d1=%f o1=%f",cn1,cd1,co1);
		sf_warning("n2=%d",cn2);
		sf_warning("Input file dimensions: ");
		sf_warning("n1=%d d1=%f o1=%f",nt,dt,ot);
		sf_warning("n2=%d d2=%f o2=%f",nh,dh,oh);
		sf_warning("n3=%d d3=%f o3=%f",nm,dm,om);
	}

	/* Read data cube */
	t=sf_floatalloc3(nt,nh,nm);
	sf_floatread(t[0][0],nh*nm*nt,in);

	/* Read cre trajectories */	
	m = sf_floatalloc2(cn1,cn2);
	sf_floatread(m[0],cn1*cn2,cremh);
	creGather = sf_floatalloc3(nt,aperture,cn2);

	mMax = om+dm*nm;
	
	for(l=0;l<nm0;l++){

		for(k=0;k<nt0;k++){
			for(i=0;i<aperture;i++){
				trac_m = (int)((double)m[(l*nt0)+k][i]/dm);

				for(j=0;j<nt;j++){
					creGather[(l*nt0)+k][i][j] = 
					(m[(l*nt0)+k][i] <= mMax)?
					t[trac_m][i][j] : 0.;
				}
			}
		}/* loop over t0s */
	}/* loop over m0s */

	/* eixo = sf_maxa(n,o,d)*/
	ax = sf_maxa(nt, ot, dt);
	ay = sf_maxa(aperture, oh, dh);
	az = sf_maxa(cn2, 0, 1);

	/* sf_oaxa(arquivo, eixo, índice do eixo) */
	sf_oaxa(out,ax,1);
	sf_oaxa(out,ay,2);
	sf_oaxa(out,az,3);
	sf_floatwrite(creGather[0][0],cn2*aperture*nt,out);

	/* eixo do vetor m */
	am1 = sf_maxa(cn1,oh,dh);
	am2 = sf_maxa(1,0,1);
	sf_oaxa(out_m,am1,1);
	sf_oaxa(out_m,az,2);
	sf_oaxa(out_m,am2,3);
	sf_floatwrite(m[0],cn1*cn2,out_m);

	exit(0);
}
