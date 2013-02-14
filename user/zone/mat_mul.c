/*
 *  mat_mul.c
 *  3D tools
 *
 *  Created by Yanadet Sripanich on 2/8/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "mat_mul.h"
#include <stdio.h>

void mat_mul(float **m1 /*input 2x2 matrix*/,float **m2 /*input 2x2 matrix*/, float **output)
/*<matrix multiplication and output>*/

{
	
	output[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0];
	output[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1];
	output[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0];
	output[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1];
	
}

