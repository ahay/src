/* Very Fast Simulated Reannealing */
/*
  Copyright (C) 2008 University of Texas at Austin
  Copyright (C) 1989-1993 Lester Ingber
  
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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vfsr_defs.h"

#ifndef _VFSR_H_
#define _VFSR_H_

double vfsr (double (*user_cost_function) (),
             double (*user_random_generator) (),
             int number_parameters, int *parameter_type,
             double *parameter_initial_final,
             double *parameter_minimum, double *parameter_maximum,
             double *tangents, double *curvature,
             int *exit_status,
             VFSR_DEFINES *OPTIONS,
             void *user_data);

/* use built-in random numbers generator */
double vfsr_std_rng (double (*user_cost_function) (),
                     int number_parameters, int *parameter_type,
                     double *parameter_initial_final,
                     double *parameter_minimum, double *parameter_maximum,
                     double *tangents, double *curvature,
                     int *exit_status,
                     VFSR_DEFINES *OPTIONS,
                     void *user_data);

#endif /* _VFSR_H_ */
