// 
// time-domain acoustic FD modeling
// 

#include <valarray>
#include <iostream>
#include <rsf.hh>
#include "cub.hh"
#include "vai.hh"
using namespace std;

int main(int argc, char* argv[])
{
    // Laplacian coefficients
    float c0=-30./12.;
    float c1=+16./12.;
    float c2=- 1./12.;

    // init RSF
    sf_init(argc,argv);

    bool verb; // vebose flag
    if(! sf_getbool("verb",&verb)) verb=0;

    CUB Fw("in", "i"); Fw.headin(); //Fw.report();
    CUB Fv("vel","i"); Fv.headin(); //Fv.report();
    CUB Fr("ref","i"); Fr.headin(); //Fr.report();
    
    // cube axes
    axis at,az,ax;
    Fw.getax(0,&at);
    Fv.getax(0,&az);
    Fv.getax(1,&ax);

    // setup output header
    CUB Fo("out","o");
    Fo.setup(3,Fv.esize());
    Fo.putax(0,&az);
    Fo.putax(1,&ax);
    Fo.putax(2,&at); //Fo.report();
    Fo.headou();

    float idx,idz,dt2;
    dt2 =    at.d*at.d ;
    idz = 1/(az.d*az.d);
    idx = 1/(ax.d*ax.d);

    // read wavelet, velocity and reflectivity
    valarray<float> ww( at.n      ); ww=0;
    valarray<float> vv( az.n*ax.n ); vv=0;
    valarray<float> rr( az.n*ax.n ); rr=0;
   
    Fw >> ww;
    Fr >> rr;
    Fv >> vv;

    // allocate temporary arrays
    valarray<float> um(az.n*ax.n); um=0;
    valarray<float> uo(az.n*ax.n); uo=0;
    valarray<float> up(az.n*ax.n); up=0;
    valarray<float> ud(az.n*ax.n); ud=0;

    // init ValArray Index counter
    VAI k(az.n,ax.n);

    // 
    // MAIN LOOP
    // 
    if(verb) cerr << endl;
    for (int it=0; it<at.n; it++) {
	if(verb) cerr << "\b\b\b\b\b" << it;

	// 4th order laplacian
	for (int iz=2; iz<az.n-2; iz++) {
	    for (int ix=2; ix<ax.n-2; ix++) {
		ud[k(iz,ix)] = 
		    c0* uo[ k(iz  ,ix  )] * (idx+idz) +
		    c1*(uo[ k(iz  ,ix-1)]+uo[ k(iz  ,ix+1)]) * idx + 
		    c1*(uo[ k(iz-1,ix  )]+uo[ k(iz+1,ix  )]) * idz + 
		    c2*(uo[ k(iz  ,ix-2)]+uo[ k(iz  ,ix+2)]) * idx + 
		    c2*(uo[ k(iz-2,ix  )]+uo[ k(iz+2,ix  )]) * idz;
	    }
	}

	ud *= vv*vv;
	
	// inject wavelet
	ud -= ww[it] * rr;
	
	// time step
	up=(float)2 * uo - um + ud * dt2;
	um =   uo;
	uo =   up;
	
	// write wavefield to output output
	Fo << uo;
    }
    if(verb) cerr << endl;

    exit(0);
}

