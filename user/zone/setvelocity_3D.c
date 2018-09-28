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

#include "ml_traveltime_vgradient_3D.h"
#include "ml_traveltime_vconstant_3D.h"
#include "ml_traveltime_vti_3D.h"
#include "general_traveltime_3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void setfunc(int vstatus, func3 *f) 
/*<Set the functions(vgradient or vconstant)>*/

{
	
	if (vstatus == 0) {
		
		f->T_k = T0_k;
		f->T_k_k_1 = T0_k_k_1;
		f->T_k_k_2 = T0_k_k_2;
		f->T_k_k1_1 = T0_k_k1_1;
		f->T_k_k1_2 = T0_k_k1_2;
		f->T_k_k_k_1 = T0_k_k_k_1;
		f->T_k_k_k_2 = T0_k_k_k_2;
		f->T_k_k_k_12 = T0_k_k_k_12;
		f->T_k_k1_k1_1 = T0_k_k1_k1_1;
		f->T_k_k1_k1_2 = T0_k_k1_k1_2;
		f->T_k_k1_k1_12 = T0_k_k1_k1_12;
		f->T_k_k_k1_1 = T0_k_k_k1_1;
		f->T_k_k_k1_2 = T0_k_k_k1_2;
		f->T_k_k_k1_12 = T0_k_k_k1_12;
		f->T_k_k_k1_21 = T0_k_k_k1_21;
		f->T_k_zk = T0_k_zk;
		f->T_k_zk1 = T0_k_zk1;
		f->T_k_zk_zk = T0_k_zk_zk;
		f->T_k_zk1_zk1 = T0_k_zk1_zk1;
		f->T_k_zk_zk1 = T0_k_zk_zk1;
		f->T_k_k_zk_1 = T0_k_k_zk_1;
		f->T_k_k_zk_2 = T0_k_k_zk_2;
		f->T_k_k1_zk1_1 = T0_k_k1_zk1_1;
		f->T_k_k1_zk1_2 = T0_k_k1_zk1_2;
		f->T_k_k_zk1_1 = T0_k_k_zk1_1;
		f->T_k_k_zk1_2 = T0_k_k_zk1_2;
		f->T_k_k1_zk_1 = T0_k_k1_zk_1;
		f->T_k_k1_zk_2 = T0_k_k1_zk_2;
		
	}
	else if (vstatus == 1) {
		
		f->T_k = T1_k;
		f->T_k_k_1 = T1_k_k_1;
		f->T_k_k_2 = T1_k_k_2;
		f->T_k_k1_1 = T1_k_k1_1;
		f->T_k_k1_2 = T1_k_k1_2;
		f->T_k_k_k_1 = T1_k_k_k_1;
		f->T_k_k_k_2 = T1_k_k_k_2;
		f->T_k_k_k_12 = T1_k_k_k_12;
		f->T_k_k1_k1_1 = T1_k_k1_k1_1;
		f->T_k_k1_k1_2 = T1_k_k1_k1_2;
		f->T_k_k1_k1_12 = T1_k_k1_k1_12;
		f->T_k_k_k1_1 = T1_k_k_k1_1;
		f->T_k_k_k1_2 = T1_k_k_k1_2;
		f->T_k_k_k1_12 = T1_k_k_k1_12;
		f->T_k_k_k1_21 = T1_k_k_k1_21;
		f->T_k_zk = T1_k_zk;
		f->T_k_zk1 = T1_k_zk1;
		f->T_k_zk_zk = T1_k_zk_zk;
		f->T_k_zk1_zk1 = T1_k_zk1_zk1;
		f->T_k_zk_zk1 = T1_k_zk_zk1;
		f->T_k_k_zk_1 = T1_k_k_zk_1;
		f->T_k_k_zk_2 = T1_k_k_zk_2;
		f->T_k_k1_zk1_1 = T1_k_k1_zk1_1;
		f->T_k_k1_zk1_2 = T1_k_k1_zk1_2;
		f->T_k_k_zk1_1 = T1_k_k_zk1_1;
		f->T_k_k_zk1_2 = T1_k_k_zk1_2;
		f->T_k_k1_zk_1 = T1_k_k1_zk_1;
		f->T_k_k1_zk_2 = T1_k_k1_zk_2;	
	}
	else if (vstatus == 2) {
		
		f->T_k = T2_k;
		f->T_k_k_1 = T2_k_k_1;
		f->T_k_k_2 = T2_k_k_2;
		f->T_k_k1_1 = T2_k_k1_1;
		f->T_k_k1_2 = T2_k_k1_2;
		f->T_k_k_k_1 = T2_k_k_k_1;
		f->T_k_k_k_2 = T2_k_k_k_2;
		f->T_k_k_k_12 = T2_k_k_k_12;
		f->T_k_k1_k1_1 = T2_k_k1_k1_1;
		f->T_k_k1_k1_2 = T2_k_k1_k1_2;
		f->T_k_k1_k1_12 = T2_k_k1_k1_12;
		f->T_k_k_k1_1 = T2_k_k_k1_1;
		f->T_k_k_k1_2 = T2_k_k_k1_2;
		f->T_k_k_k1_12 = T2_k_k_k1_12;
		f->T_k_k_k1_21 = T2_k_k_k1_21;
		f->T_k_zk = T2_k_zk;
		f->T_k_zk1 = T2_k_zk1;
		f->T_k_zk_zk = T2_k_zk_zk;
		f->T_k_zk1_zk1 = T2_k_zk1_zk1;
		f->T_k_zk_zk1 = T2_k_zk_zk1;
		f->T_k_k_zk_1 = T2_k_k_zk_1;
		f->T_k_k_zk_2 = T2_k_k_zk_2;
		f->T_k_k1_zk1_1 = T2_k_k1_zk1_1;
		f->T_k_k1_zk1_2 = T2_k_k1_zk1_2;
		f->T_k_k_zk1_1 = T2_k_k_zk1_1;
		f->T_k_k_zk1_2 = T2_k_k_zk1_2;
		f->T_k_k1_zk_1 = T2_k_k1_zk_1;
		f->T_k_k1_zk_2 = T2_k_k1_zk_2;	
	}	
	
}
