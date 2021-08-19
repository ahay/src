/* Common Reflection Element (CRE) stacking

Programmer: Rodolfo A. C. Neves (Dirack) 06/10/2019

Email:  rodolfo_profissional@hotmail.com

License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

 */

#include <rsf.h>

int main(int argc, char* argv[])
{

	float om0; // CMP axis origin
	float dm0; // CMP sampling
	int nm0; // Number of CMP's
	float oh; // Offset axis origin
	float dh; // Offset sampling
	int nh; // Number of Offsets
	int nt; // Number of time samples
	float ot; // Time axis origin
	float dt; // Time sampling
	int nt0; // Number of t0s in stacked section
	float ot0; // t0s axis origin
	float dt0; // t0s sampling
	float om0t; // CMP axis origin
	float dm0t; // CMP sampling
	int nm0t; // Number of CMP's
	float oht; // Offset axis origin
	float dht; // Offset sampling
	int nht; // Number of Offsets
	int nt0t; // Number of samples in CRE stacking curve 
	float ot0t; // Axis origin of CRE stacking curve
	float dt0t; // Sampling in CRE stacking curve
	bool verb; // Key to turn On/Off active mode
	float ***creTimeCurve; // CRE traveltime curves
	float ****creGatherCube; // CRE Data cube A(m0,t0,h,t)
	float **stackedSection; // CRE stacked section
	int it0, im0, ih, tetai; // loop counter and indexes
	float sumAmplitudes; // Sum of amplitudes in the stacking
	int aperture; // Number of offsets to stack

	/* RSF files I/O */  
	sf_file in, timeCurves, out;

	/* RSF files axis */
	sf_axis ax,ay,az;

	sf_init(argc,argv); 

	in = sf_input("in");
	timeCurves = sf_input("timeCurves");
	out = sf_output("out");

	/* Read cre gather cube geometry */
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
	if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&oh)) sf_error("No o2= in input");
	if (!sf_histint(in,"n3",&nt0)) sf_error("No n3= in input");
	if (!sf_histfloat(in,"d3",&dt0)) sf_error("No d3= in input");
	if (!sf_histfloat(in,"o3",&ot0)) sf_error("No o3= in input");
	if (!sf_histint(in,"n4",&nm0)) sf_error("No n4= in input");
	if (!sf_histfloat(in,"d4",&dm0)) sf_error("No d4= in input");
	if (!sf_histfloat(in,"o4",&om0)) sf_error("No o4= in input");

	/* Read cre time curves geometry */
	if (!sf_histint(timeCurves,"n1",&nht)) sf_error("No n1= in timeCurves input");
	if (!sf_histfloat(timeCurves,"d1",&dht)) sf_error("No d1= in timeCurves input");
	if (!sf_histfloat(timeCurves,"o1",&oht)) sf_error("No o1= in timeCurves input");
	if (!sf_histint(timeCurves,"n2",&nt0t)) sf_error("No n2= in timeCurves input");
	if (!sf_histfloat(timeCurves,"d2",&dt0t)) sf_error("No d2= timeCurves in input");
	if (!sf_histfloat(timeCurves,"o2",&ot0t)) sf_error("No o2= timeCurves in input");
	if (!sf_histint(timeCurves,"n3",&nm0t)) sf_error("No n3= in timeCurves input");
	if (!sf_histfloat(timeCurves,"d3",&dm0t)) sf_error("No d3= in timeCurves input");
	if (!sf_histfloat(timeCurves,"o3",&om0t)) sf_error("No o3= in timeCurves input");

	if(! sf_getbool("verb",&verb)) verb=0;
	/* 1: active mode; 0: quiet mode */

	if(!sf_getint("aperture",&aperture)) aperture=1;
	/* Stacking aperture, number of offsets to stack */

	if(aperture > nh){
		sf_error("The aperture can't be > n2\nAperture=%i n2=%i",aperture,nh);
	}

	if (verb) {

		sf_warning("Active mode on!!!");
		sf_warning("Input file parameters: ");
		sf_warning("n1=%i d1=%f o1=%f",nt,dt,ot);
		sf_warning("n2=%i d2=%f o2=%f",nh,dh,oh);
		sf_warning("n3=%i d3=%f o3=%f",nt0,dt0,ot0);
		sf_warning("n4=%i d4=%f o4=%f",nm0,dm0,om0);

		sf_warning("Time curves file parameters: ");
		sf_warning("n1=%i d1=%f o1=%f",nht,dht,oht);
		sf_warning("n2=%i d2=%f o2=%f",nt0t,dt0t,ot0t);
		sf_warning("n3=%i d3=%f o3=%f",nm0t,dm0t,om0t);
	}

	creTimeCurve = sf_floatalloc3(nht,nt0t,nm0t);
	sf_floatread(creTimeCurve[0][0],nm0t*nt0t*nht,timeCurves);
	creGatherCube = sf_floatalloc4(nt,nh,nt0,nm0);
	sf_floatread(creGatherCube[0][0][0],nm0*nt0*nh*nt,in);
	stackedSection = sf_floatalloc2(nt0,nm0);

	/* CRE STACKING */	
	for (im0=0; im0 < nm0; im0++){
			
		for(it0=0; it0 < nt0; it0++){

			sumAmplitudes = 0;

			for(ih=0; ih < aperture; ih++){

				tetai = (int) ((double)creTimeCurve[im0][it0][ih]/dt);
				sumAmplitudes += creGatherCube[im0][it0][ih][tetai];
				
			} /* loop over h*/

			stackedSection[im0][it0] = sumAmplitudes;

		}/* loop over t0 */

	}/* loop over m0 */

	/* axis = sf_maxa(n,o,d)*/
	ax = sf_maxa(nt0, ot0, dt0);
	ay = sf_maxa(nm0, om0, dm0);
	az = sf_maxa(1, 0, 1);

	/* sf_oaxa(file, axis, axis index) */
	sf_oaxa(out,ax,1);
	sf_oaxa(out,ay,2);
	sf_oaxa(out,az,3);
	sf_oaxa(out,az,4);
	sf_floatwrite(stackedSection[0],nt0*nm0,out);

	exit(0);
}
