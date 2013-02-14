/*
 *  setvelocity_3D.c
 *  
 *
 *  Created by Yanadet Sripanich on 2/6/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "ml_traveltime_vgradient_3D.h"
#include "ml_traveltime_vconstant_3D.h"
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
	else {
		
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
		f->T_k_k1_zk_2 = T1_k_k1_zk_2;	}	
	
}