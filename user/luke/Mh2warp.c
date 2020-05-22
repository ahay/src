/* Distance From Midpoint-squared warping. */
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

void transpose_array(float* array, float* transp, int n2, int n1){
	/* program for transposing 2d array */
	for ( int i = 0 ; i < n2 ; i++ ){
		for ( int j = 0 ; j < n1 ; j++ ){
			transp[j*n2 + i ] = array[i*n1 + j];
		}
	}
	return;
}

void get_column( float* array, float* column, int i, int n )
/* grab ith column from an array with height n */	
{
	for ( int j =  0 ; j < n ; j++ ){
		column [ j ] = array [ j + i*n ];
	}
	return;
}

void put_column( float* array, float* column, int i, int n )
/* put ith column into an array with height n */	
{
	for ( int j =  0 ; j < n ; j++ ){
		array [ j + i*n ] = column [ j ] ;
	}
	return;
}



int main(int argc, char* argv[])
{
    sf_map4 mo;
    int n1, in_n2, out_n2;
    float o1, d1, in_o2, in_d2, eps, h;
	float out_o2, out_d2, mid;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

	/* get sampling */
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
		
	if (!sf_histint(in,"n2",&in_n2)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&in_d2)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&in_o2)) sf_error("No o2= in input");

	/* input width */
	float width = in_d2*(float)(in_n2-1);
	
	/* last position on second axis */
	float in_f2 = in_o2+width;

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */	
	
	if (!sf_getint("pad",&out_n2)) out_n2=in_n2; 
	/* output time samples */

    if (!sf_getfloat("mid",&mid)) mid=(float)(in_n2-1)*in_d2/2+in_o2;
    /* center for midpoint streching */	
	if (mid < in_o2 || mid > in_f2){
		mid = in_o2+width/2;
		sf_warning("mid should be within second axis range, changing to mid=%g",mid);
	}
	/* put input sampling in output file */
	sf_putint(out,"n2_h2warp",in_n2);
	sf_putfloat(out,"d2_h2warp",in_d2);
	sf_putint(out,"o2_h2warp",in_o2);
	sf_putfloat(out,"mid_h2warp",mid);
	/* allocate array for input data */
	float* arrIn = sf_floatalloc(in_n2*n1);
	/* read data */
	sf_floatread(arrIn,in_n2*n1,in);
	/* allocate array for transposed input data */
	float* arrInT = sf_floatalloc(in_n2*n1);
	/* transpose */
	transpose_array(arrIn, arrInT, in_n2, n1);
	/* free unneded array */
	free ( arrIn);
	/* how much dist on the left side ?*/
	float left = mid - in_o2;
	/* how much dist on the right side? */
	float right = in_f2 - mid;
	/* output origin */
	out_o2 = -left*left;
	/* output sampling */
	out_d2 = (right*right + left*left)/(out_n2-1);
	/* set sampling of output file */
	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);
	sf_putint(out,"n2",out_n2);
	sf_putfloat(out,"d2",out_d2);
	sf_putfloat(out,"o2",out_o2);
   
	/* allocate intermediate arrays for stretching operations */
    float* trace = sf_floatalloc(in_n2);
    float* h2 = sf_floatalloc(out_n2);
    float* trace2 = sf_floatalloc(out_n2);
	/* initialize the stretch */
    mo = sf_stretch4_init (in_n2, in_o2, in_d2, out_n2, eps);
	
	/* get the streaching coordinates */
    for (int i =0; i  < out_n2; i++) {
		h = out_o2 + i*out_d2;
		if (h<0){
			h2[ i] = mid - sqrtf(-1*h);
		}else{
			h2[ i] = mid + sqrtf(h);
		}
    }
	/* define the stretch */
    sf_stretch4_define (mo,h2);
	/* allocate output array in transposed coords */
	float* arrOutT = sf_floatalloc(out_n2*n1);
	
	/* loop through traces of transposed array */
    for (int i = 0; i < n1; i++) {
		/* get this column */
		get_column( arrInT , trace , i, in_n2 );
		for (int j = 0 ; j < in_n2 ; j++){
		}
		/* warp the column */
		sf_stretch4_invert (false,mo,trace2,trace);
		for (int j = 0 ; j < in_n2 ; j++){
		}
		/* put it in the transposed output array */
		put_column( arrOutT, trace2, i, out_n2);
    }
	/* allocate an array for transposing */
	float* arrOut = sf_floatalloc(n1*out_n2);
	/* transpose arrOut so it has the correct axes */
	transpose_array(arrOutT, arrOut, n1, out_n2);
	/* write the warped array to disk */
	sf_floatwrite (arrOut,out_n2*n1,out);
	/* deallocate unneded arrays */
	free ( arrInT );
	free ( arrOut );
	free ( arrOutT);
	free ( h2     );
	free ( trace  );
	free ( trace2 );
    exit(0);
}
