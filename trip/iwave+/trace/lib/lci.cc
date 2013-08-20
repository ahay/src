#include "lci.h"

#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

void lci(int nx,       
	 int ny,       
	 float dy,     
	 float oy,     
	 float * yx,   
	 float * fy,   
	 float * fx    
	 ) {
  
  float denom2 = 1.0/(2.0e+00*dy*dy*dy);
  float denom6 = 1.0/(6.0e+00*dy*dy*dy);
  float maxy = oy+ny*dy;
  
  float y  = 0.0e+00;
  float ym1= 0.0e+00;
  float y0 = 0.0e+00;
  float y1 = 0.0e+00;
  float y2 = 0.0e+00;

  float p0 = 0.0e+00;
  float p1 = 0.0e+00;

  double oops;
  double duh;

  int i;
  int j = 0;
  int k = 0;

  for (i=0; i<nx; i++) {
    
    /* cannot use in c90 - requires c99
    y = fminf(maxy,fmaxf(yx[i],oy));
    */
    oops=oy;
    duh=yx[i];
    duh=MAX(duh,oops);
    oops=maxy;
    y=MIN(oops,duh);

    j = iwave_min(ny-2,iwave_max((int)((y-oy)/dy),0));
    k = iwave_min(iwave_max(0,j-1),ny-4);

    ym1 = y - k*dy;
    y0  = ym1 - dy;
    y1  = y0  - dy;
    y2  = y1  - dy;
    
    p0 = denom6*y0*y1;
    p1 = denom2*ym1*y2;

    fx[i] =
      - p0*y2*fy[k] 
      + p1*y1*fy[k+1] 
      - p1*y0*fy[k+2] 
      + p0*ym1*fy[k+3];
  }
}

