// 2-D finite-difference born modeling
#include <iostream>
#include <rsf.hh>
#include <valarray>
#include <cstddef>

using namespace std;

// Index class
class Index{
	private:
		int n[2];
	public:
		// constructor
		Index(int m0, int m1) {n[0]=m0; n[1]=m1;}
		// coordinate transformation (2D->1D)
		int operator() (int, int);
};

int Index::operator() (int i0, int i1)
{
	int ii;
	ii=i1*n[0]+i0;
	return ii;
}

static void deriv2(valarray<float> &ww, int nt)
/* second time derivative */
{
	int it;
	valarray<float> temp(nt);
	temp[0]=ww[0];
	for(it=1; it<nt; it++)
		temp[it]=ww[it]-ww[it-1];
	for(it=0; it<nt-1; it++)
		ww[it]=temp[it+1]-temp[it];
	ww[nt-1]=temp[nt-1]-temp[nt-2];
}

int main(int argc, char* argv[])
{
	// initialize Madagascar
	sf_init(argc,argv);

	// input files
	iRSF Fr, Fw("wav"), Fv("v"), Fdetv("detv");
	
	// command-line input
	iRSF par(0);

	// output file
	oRSF Fd, Fo("snapshot");

	// forward modeling or born modeling
	bool born;
	par.get("born", born, true);

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

	int nr, r0, rz, ft, jt;
	// trace number of shot record
	par.get("nr", nr, n2);
	// starting position of shot record
	par.get("r0", r0, 0);
	// depth of shot record
	par.get("rz", rz, 0);
	// first recorded time
	par.get("ft", ft, 0);
	// time interval
	par.get("jt", jt, 1);

	// set the dimension of output data file
	Fd.put("n1", nt);
	Fd.put("d1", dt);
	float zero=0.;
	Fd.put("o1", zero);
	Fd.put("label1", "Time");
	Fd.put("unit1", "s");

	Fd.put("n2", nr);
	Fd.put("d2", d2);
	Fd.put("o2", r0*d2);

	// set the dimension of output wavefield file
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

	// read wavelet, velocity, source position & reflectivity
	valarray<float> ww(nt); ww=0; Fw >> ww;
	valarray<float> vv(n12); vv=0; Fv >> vv;
	valarray<float> rr(n12); rr=0; Fr >> rr;
	valarray<float> ref(n12); ref=0; Fdetv >> ref;

	// allocate wavefield and data arrays
	size_t large=n12*nt;
	valarray<float> wave(large);
	valarray<float> dd(nt*nr);
	valarray<float> vv2(n12);

	// allocate temporary arrays 
	valarray<float> u0(n12); u0=0;
	valarray<float> u1(n12); u1=0;
	valarray<float> u2(n12); u2=0;
	valarray<float> ud(n12); ud=0;
	vv2 = vv*vv*dt2;

	// transform 2D [n1, n2] to 1D [n1*i2+i1]
	Index k(n1,n2);

	// second time derivative
	if(born) deriv2(ww, nt);

	// Time loop 1 (U_0)
	for (int it=0; it<nt; it++){
		sf_warning("Loop 1, it=%d;", it);

		// wavefield storage
		wave[slice((size_t)it*n12, n12, 1)]=u1;

		// 4th order laplacian
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
		ud *= vv2;

		//inject wavelet
		ud += ww[it]*rr;

		//time step
		u2 = (float)2*u1-u0+ud;
		u0=u1;
		u1=u2;
	} //end of it

	// if forward modeling, solve wave equation once
	if(!born){
		// write wavefield to output
		for(int it=0; it<nt; it++){
			dd[slice(it,nr,nt)]=
				wave[slice((size_t)it*n12+r0*n1+rz,nr,n1)];

			if(it >=ft && 0 == (it-ft)%jt){
				u1=wave[slice((size_t)it*n12, n12, 1)];
				Fo << u1;
			}
		}

		Fd << dd;

		// if forward modeling, only solve wave equation once
		exit(0);
	}

	// second initialization
	u0=0.; u1=0.; u2=0.; ud=0.;

	// Time loop 2 (det U)
	for (int it=0; it<nt; it++) {
		sf_warning("Loop 2, it=%d;", it);

		// write shot record to output
		dd[slice(it,nr,nt)]=
			u1[slice(r0*n1+rz,nr,n1)];

		// write wavefield to output
		if (it >= ft && 0 == (it-ft)%jt) {
			Fo << u1;
		}

		// 4th order laplacian
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

		// scale by velocity
		ud *= vv2;

		// inject source term
		u2=wave[slice((size_t)it*n12,n12,1)];
		ud += (float)2*u2*ref/vv;

		//time step
		u2 = (float)2*u1-u0+ud;
		u0=u1;
		u1=u2;
	} //end of it

	Fd << dd;

	exit(0);
}
