/* Additional defines for Very Fast Simulated Reannealing */
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

#ifndef _VFSR_DEFS_H_
#define _VFSR_DEFS_H_

#define TRUE        1
#define FALSE       0

#define INTEGER_TYPE     1
#define REAL_TYPE        0

#define NORMAL_EXIT                     0
#define P_TEMP_TOO_SMALL                1
#define C_TEMP_TOO_SMALL                2
#define COST_REPEATING                  3
#define TOO_MANY_INVALID_STATES         4

typedef struct {
    double COST_PRECISION;
    int USER_INITIAL_PARAMETERS;
    double ACCEPTED_TO_GENERATED_RATIO;
    int LIMIT_ACCEPTANCES;
    double TEMPERATURE_RATIO_SCALE;
    double TEMPERATURE_ANNEAL_SCALE;
    double COST_PARAMETER_SCALE;
    int TESTING_FREQUENCY_MODULUS;
    int MAXIMUM_REANNEAL_INDEX;
    double REANNEAL_RESCALE;
    double INITIAL_PARAMETER_TEMPERATURE;
    int USER_INITIAL_PARAMETERS_TEMPS;
    int USER_INITIAL_COST_TEMP;
    int NUMBER_COST_SAMPLES;
    int MAXIMUM_COST_REPEAT;
    double DELTA_X;
    int INCLUDE_INTEGER_PARAMETERS;
    int ACTIVATE_REANNEAL;
    int LIMIT_INVALID_GENERATED_STATES;
} VFSR_DEFINES;

#endif /* _VFSR_DEFS_H_ */
