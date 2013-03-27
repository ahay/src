#include <rsf.h>

//Program Computes the distance a point is from a given elipse
//Nearest point to an elipse for items in the first quadrant
int main (int argc, char* argv[]) 
{
//declare variables
   float y0, y1, x0, x1, e0, e1, esqr1, esqr0, ey0, ey1, t, t0, t1, r0, r1, f, d0, d1, distance, e0y0, x0de0, x0de0sqr, denom0;



	if (y1 > 0){
		
		if (y0 > 0){
			//Bisect to compute root for parameter t
			esqr0 = e0*e0;
                        esqr1 = e1*e1;
  			ey0 = e0*y0;
  			ey1 = e1*y1;
                        
                        t0 = -esqr1 + ey1;
                        t1 = -esqr1 + sqrt(ey0*ey0+ey1*ey1);
                        t = t0;
		const int imax = 2 * FLT_MAX_EXP;
			for (int i=0; i < imax; i++){
			t = (0.5)*(t0+t1);
                        	if (t == t0 || t == t1){
				break;
				}			
			r0 = ey0/(t + esqr0);
			r1 = ey1/(t + esqr1);

			f = r0*r0 + r1*r1 -1.;
				if (f > 0.){
				t0=t;
				}else if (f < 0.){
				t1=t;
				}else{
				break;
				}
			}
			x0 = esqr0*y0/(t+esqr0);
			x1 = esqr1*y1/(t+esqr1);
			d0 = x0 - y0;
			d1 = x1 - y1;

			distance = sqrt(d0*d0 + d1*d1);
		} else { //y0 == 0
			x0=0.;
			x1=e1;
			distance = fabs(y1-e1);
		}
	}else{ //y1==0
		denom0 = e0*e0 - e1*e1;
		e0y0 = e0*y0;
		if (e0y0 < denom0){
		// y0 is inside the subinterval
		x0de0 = e0y0/denom0;
		x0de0sqr = x0de0*x0de0;
		x0 = e0*x0de0;
		x1 = e1*sqrt(fabs(1. - x0de0sqr));
		d0 = x0 - y0;
		distance = sqrt(d0*d0 + x1*x1);
		}else{
		// y0 is outside the subinterval.  The closest ellipse point has
		//x1 ==0 and is on the boundary interval
		x0 = e0;
		x1 = 0.;
		distance = fabs(y0 - e0);
	}
	return distance;
}

}





