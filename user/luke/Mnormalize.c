/* Normalization of data batches up to a given axis.
   Can either normalize to be between [0,1] with type=m
   or have mean=0, std dev = 1 with type=s */
/*
  Copyright (C) 2020 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#include <math.h>

float get_minimum(float* array, int n){
    /* returns the minimum value of an array */
    int i;
    float min;
    min = array[0];
    for (i = 0 ; i < n ; i++){
        /* make sure value is not NaN */
        if (array[i] != array[i]) continue;
	      if (array[i] < min) min=array[i];
    }
    return min;
}

float get_maximum(float* array, int n){
    /* returns the maximum value of an array */
    int i;
    float max;
    max = array[0];
    for (i = 0 ; i < n ; i++){
        /* make sure value is not NaN */
        if (array[i] != array[i]) continue;
	      if (array[i] > max) max=array[i];
    }
    return max;
}

float get_mean(float* array, int n){
    /* returns the mean of an array */
    int i;
    double mean;
    mean = 0.0;
    for (i = 0 ; i < n ; i++){
        /* make sure value is not NaN */
        if (array[i] != array[i]) continue;
        mean += (double)array[i];
    }
    return (float)(mean/(double)n);
}

float get_std(float* array, float mean, int n){
    /* returns the standard deviation of an array given a mean value */
    int i;
    double std, dmean;
    std = 0.0;
    dmean = (double)mean; 
    for (i = 0 ; i < n ; i++){
         if (array[i] != array[i]) continue;
         std += ((double)array[i]-dmean)*((double)array[i]-dmean);
    }
    return (float) sqrt(std/(double)n);
}

int main(int argc, char* argv[])
{
    int i, i1, i2, axis, ndim;
    off_t n1, n2, n[SF_MAX_DIM];
    float *trace=NULL;
    sf_file in=NULL, out=NULL; /* Input and output files */
    float mean, std, min, max;
    float bias, scale;
    char *type;

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in))
	sf_error("Need float input");

    if (!sf_getint("axis",&axis)) axis = 1;
    /* normalize so axes groups up to this dimension have mean 0, std deviation 1 */

    if (NULL == (type = sf_getstring("type"))) type = "s";
    /* 'm' means data are biased and scaled to be between 0 and 1 
       's' means data are biased and scaled to have mean=0, std dev = 1  
    */

    /* how many dimensions?*/
    ndim = sf_largefiledims (in, n);

    n1=1;
    n2=1;
    /* assign dimensions lower than axis to n1, larger to n2 */
    for (i=0; i < ndim; i++) {
	if (i < axis) n1 *= n[i];
	else          n2 *= n[i];
    }

    /* allocate floating point array */
    trace = sf_floatalloc (n1);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	/*read a trace */
	sf_floatread(trace,n1,in);
	/* get min and max of trace */
        min = get_minimum(trace, n1);
        max = get_maximum(trace, n1);
        switch (type[0]){
            case 'm':
               /* 0,1 normalization */
	       /* get trace min */
               min = get_minimum(trace, n1);
	       /* get trace max */
               max = get_maximum(trace, n1);
	       /* write to bias and std */
               bias = min;
               scale = max-min;
               break;
            case 's':
            default:
	       /* using mean zero std dev 1 normalization */
	       /*compute trace mean */
               mean = get_mean(trace, n1);
	       /* compute trace standard deviation */
	       std = get_std(trace, mean, n1);
	       /* write mean and std to bias and scale */
	       bias = mean;
	       scale = std;
	       break;
	    }
        /* address division by zero */
         if (scale == 0.0) scale = 1.0; 
	/* loop over samples */
	for (i1=0; i1 < n1; i1++) {
	    if (trace[i1] != trace[i1]) trace[i1] = 0.0;
            trace[i1] = (trace[i1] - bias)/scale;
	}

	/* write a trace */
	sf_floatwrite(trace,n1,out);
    }

    free(trace);
    sf_close();
    exit(0);
}
