// 2-D finite-difference acoustic wave propagation
#include <iostream>
#include <rsf.hh>
#include <valarray>

using namespace std;

// Index class
class Index{
	private:
		int n[2];
	public:
		// constructor
		Index(int m0, int m1) { n[0]=m0; n[1]=m1;}
		// coordinate change (2D->1D)
		int operator() (int, int);
};

int Index::operator() (int i0, int i1)
{
	int ii;
	ii=i1*n[0]+i0;
	return ii;
}

int main(int argc, char* argv[])
{
	// initialize Madagascar
	sf_init(argc,argv);

	// input files
	iRSF Fr, Fw("wav"), Fv("v");
	// command-line input
	iRSF par(0);

	// output file
	oRSF Fo;

	// Read/Write axes information
	int n1, n2, nt;
	float d1, d2, dt;
	Fr.get("n1", n1);
	Fr.get("n2", n2);
	Fr.get("d1", d1);
	Fr.get("d2", d2);
	Fw.get("n1", nt);
	Fw.get("d1", dt);

	int n12;
	n12=n1*n2;

	int ft, jt;
	// first recorded time
	par.get("ft", ft, 0);
	// time interval
	par.get("jt", jt, 1);

	// set the dimension of output file
	Fo.put("n3", (nt-ft)/jt);
	Fo.put("d3", jt*dt);
	Fo.put("o3", ft*dt);

	float dt2;
	dt2=dt*dt;

	// set laplacian coefficients
	d1=1.0/(d1*d1);
	d2=1.0/(d2*d2);

	const float c11=4.0*d1/3.0;
	const float c12=-d1/12.0;
	const float c21=4.0*d2/3.0;
	const float c22=-d2/12.0;
	const float c0=-2.0 * (c11+c12+c21+c22);

	// read wavelet, velocity & source position
	valarray<float> ww(nt); ww=0; Fw >> ww;
	valarray<float> vv(n1*n2); vv=0; Fv >> vv;
	valarray<float> rr(n1*n2); rr=0; Fr >> rr;

	// allocate temporary arrays 
	valarray<float> u0(n1*n2); u0=0;
	valarray<float> u1(n1*n2); u1=0;
	valarray<float> u2(n1*n2); u2=0;
	valarray<float> ud(n1*n2); ud=0;
	vv *= vv*dt2;

	// change 2D [n1, n2] to 1D [n1*i2+i1]
	Index k(n1,n2);

	// Time loop
	for (int it=0; it<nt; it++){

		//4th order laplacian
		for(int i1=2; i1<n1-2; i1++){
			for(int i2=2; i2<n2-2; i2++){
				ud[k(i1,i2)] = 
					c11*(u1[k(i1-1,i2)]+u1[k(i1+1,i2)]) +
					c12*(u1[k(i1-2,i2)]+u1[k(i1+2,i2)]) +
					c21*(u1[k(i1,i2-1)]+u1[k(i1,i2+1)]) +
					c22*(u1[k(i1,i2-2)]+u1[k(i1,i2+2)]) +
					c0*u1[k(i1,i2)];
			}
		}

		//scale by velocity
		ud *= vv;

		//inject wavelet
		ud += ww[it]*rr;

		//time step
		u2 = (float)2*u1-u0+ud;
		u0=u1;
		u1=u2;

		// write wavefield to output
		if(it>=ft && 0==(it-ft)%jt){
			cerr << "\b\b\b\b\b" << it;
			Fo << u1;
		}
	}
	cerr << endl;

	exit(0);
}
