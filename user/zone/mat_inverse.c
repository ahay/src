/*
 *  mat_inverse.c
 *  3D tools
 *
 *  Created by Yanadet Sripanich on 2/8/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "mat_inverse.h"
#include <stdio.h>

void mat_inverse(float **m /*input 2x2 matrix*/)
/*<Find a 2x2 matrix inverse>*/
{
	float c;
	float det;
	
	det = m[0][0]*m[1][1]-m[0][1]*m[1][0];
	
	if (det!=0) {
		m[0][0] /= det;
		m[1][1] /= det;
		m[0][1] /= -det;
		m[1][0] /= -det;
		
		c = m[0][0];
		m[0][0] = m[1][1];
		m[1][1] = c;
	}
	else {
		if (m[0][0]!=0) {
			m[0][0] = 1/m[0][0];
		}
		if (m[1][1]!=0) {
			m[1][1] = 1/m[1][1];
		}
	}

}
