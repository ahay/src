/* simple interpolator */

typedef float (*sinterp)(float *in, float x, int n);
/* generic interpolation interface */
/*^*/

static float interp0(float *in, float x, int n)
/*< nearest interpolation >*/
{
	if(x<0) return in[0];
	else if (x+0.5 >= n) return in[n-1];
	else return in[(int)(x+0.5)];
}

static float interp1(float *in, float x, int n)
/*< linear interpolation >*/
{
	int k;
	float d;

	if(x<0) return in[0];
	else if (x > n-1) return in[n-1];
	k = x;
	d = x-k;
	return ((1.0-d)*in[k] + d*in[k+1]);
}



