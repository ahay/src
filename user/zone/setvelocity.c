/*
 *  setvelocity->c
 *  
 *
 *  Created by Yanadet Sripanich on 1/17/13->
 *  Copyright 2013 __MyCompanyName__-> All rights reserved->
 *
 */
#include "ml_traveltime_vgradient.h"
#include "ml_traveltime_vconstant.h"
#include "general_traveltime.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void setfunc(int vstatus, func3 *f) 
/*<Set the functions(vgradient or vconstant)>*/

{
	
	if (vstatus == 0) {
		
		f->T_k = T0_k;
		f->T_k_k = T0_k_k;
		f->T_k_k1 = T0_k_k1;
		f->T_k_k_k = T0_k_k_k;
		f->T_k_k1_k1 = T0_k_k1_k1;
		f->T_k_k_k1 = T0_k_k_k1;
		f->T_k_zk = T0_k_zk;
		f->T_k_zk1 = T0_k_zk1;
		f->T_k_zk_zk = T0_k_zk_zk;
		f->T_k_zk1_zk1 = T0_k_zk1_zk1;
		f->T_k_zk_zk1 = T0_k_zk_zk1;
		f->T_k_k_zk = T0_k_k_zk;
		f->T_k_k1_zk1 = T0_k_k1_zk1;
		f->T_k_k_zk1 = T0_k_k_zk1;
		f->T_k_k1_zk = T0_k_k1_zk;
		
	}
	else {
		
		f->T_k = T1_k;
		f->T_k_k = T1_k_k;
		f->T_k_k1 = T1_k_k1;
		f->T_k_k_k = T1_k_k_k;
		f->T_k_k1_k1 = T1_k_k1_k1;
		f->T_k_k_k1 = T1_k_k_k1;
		f->T_k_zk = T1_k_zk;
		f->T_k_zk1 = T1_k_zk1;
		f->T_k_zk_zk = T1_k_zk_zk;
		f->T_k_zk1_zk1 = T1_k_zk1_zk1;
		f->T_k_zk_zk1 = T1_k_zk_zk1;
		f->T_k_k_zk = T1_k_k_zk;
		f->T_k_k1_zk1 = T1_k_k1_zk1;
		f->T_k_k_zk1 = T1_k_k_zk1;
		f->T_k_k1_zk = T1_k_k1_zk;
	}	
	
}