#include <math.h>

#include<rsf.h>

#include "fastmarch.h"

int main (int argc,char* argv[]) 
{
    int b1, b2, b3, n1, n2, n3, i, nshot, ndim, is,order,n123, *p;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3, slow;
    float **s, *t, *v;
    bool isvel;
    sf_file vel, time, shots;

    sf_init (argc, argv);
    vel = sf_input("in");
    time = sf_output("out");

    if (SF_FLOAT != sf_gettype(vel)) 
	sf_error("Need float input");
    if(!sf_histint(vel,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(vel,"n2",&n2)) sf_error("No n2= in input");
    if(!sf_histint(vel,"n3",&n3)) n3=1;

    if(!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in input");
    if(!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in input");
    if(!sf_histfloat(vel,"d3",&d3)) d3=d2;

    if(!sf_histfloat(vel,"o1",&o1)) o1=0.;
    if(!sf_histfloat(vel,"o2",&o2)) o2=0.;
    if(!sf_histfloat(vel,"o3",&o3)) o3=0.;

    if(!sf_getbool("vel",&isvel)) isvel=true;
    if(!sf_getint("order",&order)) order=2;

    if(!sf_getfloat("br1",&br1)) br1=d1; 
    if(!sf_getfloat("br2",&br2)) br2=d2; 
    if(!sf_getfloat("br3",&br3)) br3=d3; 
    if(!sf_getint("b1",&b1)) b1= (int) (br1/d1+0.5); 
    if(!sf_getint("b2",&b2)) b2= (int) (br2/d2+0.5); 
    if(!sf_getint("b3",&b3)) b3= (int) (br3/d3+0.5); 
    if( b1<1 ) b1=1;  
    if( b2<1 ) b2=1;  
    if( b3<1 ) b3=1;

    if(NULL != sf_getstring("shotfile")) {
	shots = sf_input("shotfile");

	if (SF_FLOAT != sf_gettype(shots)) 
	    sf_error("Need float shotfile");
	if(!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
	if(!sf_histint(shots,"n1",&ndim) || ndim != 3) 
	    sf_error("Need n1=3 in shotfile");
  
	s = sf_floatalloc2 (ndim,nshot);
	sf_read(s[0],sizeof(float),nshot*ndim,shots);
    
	sf_putint (time,"n4",nshot);
    } else {
	nshot = 1;
	ndim = 3;
    
	s = sf_floatalloc2 (ndim,nshot);     

	if(!sf_getfloat("zshot",s[0])  ) s[0][0]=0.; 
	if(!sf_getfloat("yshot",s[0]+1)) s[0][1]=o2 + 0.5*(n2-1)*d2;
	if(!sf_getfloat("xshot",s[0]+2)) s[0][2]=o3 + 0.5*(n3-1)*d3;
    
	sf_warning("Shooting from zshot=%g yshot=%g xshot=%g",
		   s[0][0],s[0][1],s[0][2]);
    }

    n123 = n1*n2*n3;

    t  = sf_floatalloc (n123);
    v  = sf_floatalloc (n123);
    p  = sf_intalloc (n123);

    sf_read(v,sizeof(float),n123,vel);
   /* transform velocity to slowness squared */
    if (isvel) {
	for(i = 0; i < n123; i++) {
	    slow = v[i];
	    v[i] = 1./(slow*slow);
	}
    } 
    
    fastmarch_init (n3,n2,n1);
  
    /* loop over shots */
    for( is = 0; is < nshot; is++) {
	fastmarch(t,v,p,
		  n3,n2,n1,
		  o3,o2,o1,
		  d3,d2,d1,
		  s[is][2],s[is][1],s[is][0], 
		  b3,b2,b1,
		  order); 
	
	sf_write (t,sizeof(float),n123,time);
    }
    
    exit (0);
}
