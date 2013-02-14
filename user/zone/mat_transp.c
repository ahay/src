/*
 *  mat_transp.c
 *  3D tools
 *
 *  Created by Yanadet Sripanich on 2/8/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include "mat_transp.h"

void mat_transp(float **m/*input 2x2 matrix*/,float **output)
/*<Find a 2x2 matrix transpose>*/
{
	
	output[0][0] = m[0][0];
	output[1][0] = m[0][1];
	output[0][1] = m[1][0];
	output[1][1] = m[1][1];
	
}
