#include <rsf.h>


float sf_traveltime(int type, float t0, float x, float v)
/*< travel time:  type =	0 line, 1 parabola, x hyperbola >*/
{
	float deltat;
	deltat= x/v;
	switch(type)
	{
	case 0:
		return (t0+deltat);
		break;
	case 1:
		return (t0+deltat*deltat);
		break;
	default:
		return sqrtf(t0*t0+deltat*deltat);
	}
}
