/*
 *  vectorsub.c
 *  Multi-Layered
 *
 *  Created by Yanadet Sripanich on 8/3/12.
 *  
 *
 */
/*
 Copyright (C) 2009 University of Texas at Austin
 
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

#include "vectorsub.h"
#include <stdio.h>
#include <math.h>
#include <rsf.h>


void vector_sub(int count1, float *inp /*y_old*/,
				int count2, float *inp_d /*d(changes)*/,
				float *outp, int count3 /*starting form which element*/)
/*<Run Newton's method for finding a good set of intersection points y>*/

{
	int n; /*number of interfaces*/	
	int i/*,j*/; /*counter*/
	
	if (count1 != count2) {
		printf("The size of the two vectors don't match");
	}
	else {
		n = count1;
	}
	
	for(i=count3;i<n;i++){
		*(outp+i) = *(inp+i) - *(inp_d+i); /*adding step*/
	}
	
	
	/*for (j=0; j<n; j++) { printf("%3.3g ", *(outp+j));} To print output if prefer*/
	
}
