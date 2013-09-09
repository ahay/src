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

/*For use in kirmod_newton */

#include <stdio.h>
#include "mat_transp.h"

void mat_transp(float **m/* Input 2x2 matrix*/,float **output)
/*<Find a 2x2 matrix transpose>*/
{
	
	output[0][0] = m[0][0];
	output[1][0] = m[0][1];
	output[0][1] = m[1][0];
	output[1][1] = m[1][1];
	
}
