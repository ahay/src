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

#include "vfsr.h"

#define ZERO ((double) 0.0)
#define ONE ((double) 1.0)
#define TWO ((double) 2.0)
#define TEN ((double) 10.0)
#define HALF ((double) 0.5)

/* Machine Options */

#if INT_LONG
#define LONG_INT long int
#else
#define LONG_INT int
#endif

/* You can define SMALL_FLOAT to better correlate to your machine's
   precision, i.e., as used in vfsr */
#ifndef SMALL_FLOAT
#define SMALL_FLOAT 1.0E-18
#endif

/* You can define your machine's maximum and minimum reals here */
#ifndef MIN_FLOAT
#define MIN_FLOAT (-1.0/SMALL_FLOAT)
#endif

/* #define MAX_FLOAT MAX_REAL */
#ifndef MAX_FLOAT
#define MAX_FLOAT (1.0/SMALL_FLOAT)
#endif

VFSR_DEFINES *OPTIONS;

/* essential MACROS */

 /* VFOR
    is a simple macro to iterate on each parameter index. */
#define VFOR(index_v) \
 for (index_v = 0; index_v < number_parameters; ++index_v)

 /* EXPONENT_CHECK
    checks that an exponent x is within a valid range and,
    if not, reduces its magnitude to fit in the range. */
#define EXPONENT_CHECK(x) \
 ((x) < min_exponent ? min_exponent : \
 ((x) > max_exponent ? max_exponent : (x)))

 /* PARAMETER_RANGE_TOO_SMALL(x)
    checks if the range of parameter x is too small to work with
    If user_cost_function changes the parameter ranges,
    this test could be (and has been) used to adaptively bypass
    some parameters, e.g., depending on constraints. */
#define PARAMETER_RANGE_TOO_SMALL(x) \
 (fabs(parameter_minimum[x] - parameter_maximum[x]) < (double) SMALL_FLOAT)

 /* INTEGER_PARAMETER(x)
    determines if the parameter is an integer type. */
#define INTEGER_PARAMETER(x) (parameter_type[x] == INTEGER_TYPE)

 /* ROW_COL_INDEX (i, j)
    converts from row i, column j to an index. */
#define ROW_COL_INDEX(i, j) (i + number_parameters * j)

 /* The State of the system, its parameters and their resulting
    function value */
typedef struct {
    double cost;
    double *parameter;
}

STATE;

 /* The 3 states that are kept track of during the annealing process */
static STATE current_generated_state, last_saved_state, best_generated_state;

 /* The array of parameter bounds */
static double *parameter_minimum, *parameter_maximum;

 /* The array of tangents (absolute value of the numerical derivatives),
    and the maximum |tangent| of the array */
static double *tangents, maximum_tangent;

 /* The parameter curvature/covariance about the best state */
static double *curvature;

 /* ratio of acceptances to generated points - determines when to
    test/reanneal */
static double accepted_to_generated_ratio;

 /* temperature parameters */
static double temperature_scale, temperature_scale_parameters;
 /* relative scalings of cost and parameters to temperature_scale */
static double temperature_scale_cost;
static double *current_parameter_temperature, *initial_parameter_temperature;
static double current_cost_temperature, initial_cost_temperature;
static double temperature_rescale_power;        /* used when applying REANNEAL_RESCALE */
static double log_new_temperature_ratio;        /* current temp = initial temp *
                                            exp(log_new_temperature_ratio) */
static double one_over_number_parameters;       /* 1/number_parameters */
static int index_exit_v;                /* information for vfsr_exit */

 /* flag to determine if curvatures should be calculated */
static int curvature_flag;

 /* counts of generated states and acceptances */
static LONG_INT *index_parameter_generations;
static LONG_INT number_generated, best_number_generated_saved;
static LONG_INT recent_number_generated, number_parameters, number_accepted;
static LONG_INT recent_number_acceptances, index_cost_acceptances;
static LONG_INT number_acceptances_saved, best_number_accepted_saved;

/* Flag indicates that the parameters generated were
   invalid according to the cost function validity criteria. */
static int valid_state_generated_flag;
LONG_INT number_invalid_generated_states, repeated_invalid_states;

/* parameter type is real or integer */
static int *parameter_type;

/* used by EXPONENT_CHECK */
static double max_exponent, min_exponent;

/* used to index repeated and recursive calls to vfsr */
static int vfsr_open = FALSE;
static int number_vfsr_open = 0;
static int recursive_vfsr_open = 0;

/* passed to user_cost_function */
static void *user_data;

 /* global parameterization of DEFINE_OPTIONS */
static double COST_PRECISION;
static int USER_INITIAL_PARAMETERS;
static double ACCEPTED_TO_GENERATED_RATIO;
static int LIMIT_ACCEPTANCES;
static double TEMPERATURE_RATIO_SCALE;
static double TEMPERATURE_ANNEAL_SCALE;
static double COST_PARAMETER_SCALE;
static int TESTING_FREQUENCY_MODULUS;
static int MAXIMUM_REANNEAL_INDEX;
static double REANNEAL_RESCALE;
static double INITIAL_PARAMETER_TEMPERATURE;
static int USER_INITIAL_PARAMETERS_TEMPS;
static int USER_INITIAL_COST_TEMP;
static int NUMBER_COST_SAMPLES;
static int MAXIMUM_COST_REPEAT;
static double DELTA_X;
static int INCLUDE_INTEGER_PARAMETERS;
static int ACTIVATE_REANNEAL;
static int LIMIT_INVALID_GENERATED_STATES;

 /* vfsr function prototypes */
static void accept_new_state(double (*user_random_generator) ());
static void generate_new_state(double (*user_random_generator) ());
static void reanneal(void);
static void cost_derivatives(double (*user_cost_function) ());
static double generate_vfsr_state(double temp,
                                  double (*user_random_generator) ());
static void vfsr_exit(double (*user_cost_function) (),
                      double *parameter_initial_final,
                      double *final_cost,
                      int *exit_status);

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: vfsr
* Parameters:
* Inputs:
*     double (*user_cost_function) ()    - user's cost function
*     double (*user_random_generator) () - user's rng
*     int g_number_parameters            - number of parameters
*     int *g_parameter_type              - type of parameter
*     double *parameter_initial_final    - initial and final parameters
*     double *g_parameter_minimum        - minimum parameter values
*     double *g_parameter_maximum        - maximum parameter values
*     double *g_tangents                 - holds 1st derivatives
*     double *g_curvature                - holds 2nd derivatives
*     int *exit_status                   - reason for exiting
*     VFSR_DEFINES *g_OPTIONS            - algorithm tuning
* Outputs:
*     double *parameter_initial_final    - holds final parameters
* Global Variables:
* Returns:
*     double final_cost                  - best cost found
* Calls:
* Description:
*        This procedure implements the full VFSR function optimization.
* Local Variables:
*        int index_cost_constraint       - index initial cost functions
*                                        - to set initial temperature
*        int temp_var_int                - integer temporary value
*        int index_v                     - iteration index
*        double index_cost_repeat        - test at MAXIMUM_COST_REPEAT
*        double sampled_cost_sum         - temporary initial costs
*        double tmp_var_db1, tmp_var_db2 - temporary doubles
* Portability:
* Other:
***********************************************************************/
double vfsr(                        /* return final_cost */
               double (*user_cost_function) (),    /* user's cost function */
               double (*user_random_generator) (), /* user's random
                                                      generator on {0,1} */
               int g_number_parameters,            /* number of parameters */
               int *g_parameter_type,              /* INTEGER_TYPE, REAL_TYPE */
               double *parameter_initial_final,    /* initial and final
                                                      parameters */
               double *g_parameter_minimum,        /* minimum parameter values */
               double *g_parameter_maximum,        /* maximum parameter values */
               double *g_tangents,         /* numerical 1st derivatives of cost */
               double *g_curvature,        /* numerical 2nd derivatives of cost */
               int *exit_status,   /* exit code */
               VFSR_DEFINES * g_OPTIONS /* [optional] */,
               void *g_user_data /* passed to user_cost_function*/)
{
    int index_cost_constraint,        /* index initial cost functions averaged
                                         to set initial temperature */
        temp_var_int,                 /* integer temporary value */
        index_v;                      /* iteration index */
    double index_cost_repeat,         /* test COST_PRECISION when this reaches
                                         MAXIMUM_COST_REPEAT */
           sampled_cost_sum,          /* temporary initial costs */
           final_cost,                /* best cost to return to user */
           tmp_var_db1, tmp_var_db2;  /* temporary doubles */

    /* initializations */

/* g_... are defined to make local information passed to vfsr() by
   the user into global variables useful for internal vfsr() routines.
   This was considered a reasonable tradeoff to passing several
   parameters between several routines. */

    number_parameters = g_number_parameters;
    parameter_type = g_parameter_type;
    parameter_minimum = g_parameter_minimum;
    parameter_maximum = g_parameter_maximum;
    tangents = g_tangents;
    curvature = g_curvature;
    OPTIONS = g_OPTIONS;
    user_data = g_user_data;

    if (OPTIONS == (VFSR_DEFINES *) NULL) {
        *exit_status = -1;
        vfsr_exit(user_cost_function,
                  parameter_initial_final,
                  &final_cost,
                  exit_status);
        return (-1);
    }
    COST_PRECISION =
        (*OPTIONS).COST_PRECISION;
    USER_INITIAL_PARAMETERS =
        (*OPTIONS).USER_INITIAL_PARAMETERS;
    ACCEPTED_TO_GENERATED_RATIO =
        (*OPTIONS).ACCEPTED_TO_GENERATED_RATIO;
    LIMIT_ACCEPTANCES =
        (*OPTIONS).LIMIT_ACCEPTANCES;
    TEMPERATURE_RATIO_SCALE =
        (*OPTIONS).TEMPERATURE_RATIO_SCALE;
    TEMPERATURE_ANNEAL_SCALE =
        (*OPTIONS).TEMPERATURE_ANNEAL_SCALE;
    COST_PARAMETER_SCALE =
        (*OPTIONS).COST_PARAMETER_SCALE;
    TESTING_FREQUENCY_MODULUS =
        (*OPTIONS).TESTING_FREQUENCY_MODULUS;
    MAXIMUM_REANNEAL_INDEX =
        (*OPTIONS).MAXIMUM_REANNEAL_INDEX;
    REANNEAL_RESCALE =
        (*OPTIONS).REANNEAL_RESCALE;
    INITIAL_PARAMETER_TEMPERATURE =
        (*OPTIONS).INITIAL_PARAMETER_TEMPERATURE;
    USER_INITIAL_PARAMETERS_TEMPS =
        (*OPTIONS).USER_INITIAL_PARAMETERS_TEMPS;
    USER_INITIAL_COST_TEMP =
        (*OPTIONS).USER_INITIAL_COST_TEMP;
    NUMBER_COST_SAMPLES =
        (*OPTIONS).NUMBER_COST_SAMPLES;
    MAXIMUM_COST_REPEAT =
        (*OPTIONS).MAXIMUM_COST_REPEAT;
    DELTA_X =
        (*OPTIONS).DELTA_X;
    INCLUDE_INTEGER_PARAMETERS =
        (*OPTIONS).INCLUDE_INTEGER_PARAMETERS;
    ACTIVATE_REANNEAL =
        (*OPTIONS).ACTIVATE_REANNEAL;
    LIMIT_INVALID_GENERATED_STATES =
        (*OPTIONS).LIMIT_INVALID_GENERATED_STATES;

    if (vfsr_open == FALSE) {
        vfsr_open = TRUE;
        ++number_vfsr_open;
    } else {
        ++recursive_vfsr_open;
    }

    /* set indices and counts to 0 */
    best_number_generated_saved =
        number_generated =
        recent_number_generated =
        recent_number_acceptances = 0;
    index_cost_acceptances =
        best_number_accepted_saved =
        number_accepted =
        number_acceptances_saved = 1;
    index_cost_repeat = 0;

    /* do not calculate curvatures initially */
    curvature_flag = FALSE;

    /* allocate storage for all parameters */
    current_generated_state.parameter =
        (double *) calloc(number_parameters, sizeof(double));
    last_saved_state.parameter =
        (double *) calloc(number_parameters, sizeof(double));
    best_generated_state.parameter =
        (double *) calloc(number_parameters, sizeof(double));
    current_parameter_temperature =
        (double *) calloc(number_parameters, sizeof(double));
    initial_parameter_temperature =
        (double *) calloc(number_parameters, sizeof(double));
    index_parameter_generations =
        (LONG_INT *) calloc(number_parameters, sizeof(LONG_INT));

    /* set the user's initial parameters to the initial system state */
    VFOR(index_v) {
        current_generated_state.parameter[index_v] =
            parameter_initial_final[index_v];
    }

    VFOR(index_v) {
        /* set all starting and current temperatures */
        if (USER_INITIAL_PARAMETERS_TEMPS) {
            current_parameter_temperature[index_v] =
                initial_parameter_temperature[index_v] =
                tangents[index_v];
        } else {
            current_parameter_temperature[index_v] =
                initial_parameter_temperature[index_v] =
                ((double) INITIAL_PARAMETER_TEMPERATURE);
        }

        /* set all states to the initial parameter values */
        best_generated_state.parameter[index_v] =
            last_saved_state.parameter[index_v] =
            current_generated_state.parameter[index_v];
    }

/* define variables for EXPONENT_CHECK */
    max_exponent = 2.302585 * log((double) MAX_FLOAT);
    min_exponent = -2.302585 * log(fabs((double) MIN_FLOAT));

    /* set the initial index of parameter generations to 1 */
    VFOR(index_v) {
        /* ignore parameters that have too small a range */
        if (PARAMETER_RANGE_TOO_SMALL(index_v))
            continue;
        index_parameter_generations[index_v] = 1;
    }

    if (USER_INITIAL_COST_TEMP) {
        initial_cost_temperature = current_cost_temperature =
            curvature[0];
    } else {
        /* calculate the average cost over samplings of the cost function */
        sampled_cost_sum = ZERO;
        for (index_cost_constraint = 0;
             index_cost_constraint < ((int) NUMBER_COST_SAMPLES);
             ++index_cost_constraint) {
            number_invalid_generated_states = 0;
            repeated_invalid_states = 0;
            do {
                ++number_invalid_generated_states;
                ++repeated_invalid_states;
                valid_state_generated_flag = TRUE;
                generate_new_state(user_random_generator);
                current_generated_state.cost =
                    (*user_cost_function)
                    (current_generated_state.parameter,
                     parameter_minimum,
                     parameter_maximum,
                     &valid_state_generated_flag,
                     OPTIONS,
                     user_data);
                if (repeated_invalid_states >
                    (int) LIMIT_INVALID_GENERATED_STATES) {
                    *exit_status = (int) TOO_MANY_INVALID_STATES;
                    vfsr_exit(user_cost_function,
                              parameter_initial_final,
                              &final_cost,
                              exit_status);
                    return (final_cost);
                }
            } while (valid_state_generated_flag == FALSE);
            --number_invalid_generated_states;
            sampled_cost_sum += fabs(current_generated_state.cost);
        }

        /* set all cost temperatures to the average cost */
        initial_cost_temperature = current_cost_temperature =
            sampled_cost_sum / ((double) NUMBER_COST_SAMPLES);
    }

    /* establish an initial current state */

    /* user's initial guess of parameters */
    VFOR(index_v) {
        current_generated_state.parameter[index_v] =
            parameter_initial_final[index_v];
    }

    if (USER_INITIAL_PARAMETERS) {        /* if using user's initial parameters */
        valid_state_generated_flag = TRUE;
        current_generated_state.cost =
            (*user_cost_function)
            (current_generated_state.parameter,
             parameter_minimum,
             parameter_maximum,
             &valid_state_generated_flag,
             OPTIONS,
             user_data);
    } else {                        /* let vfsr generate valid initial parameters */
        repeated_invalid_states = 0;
        do {
            ++number_invalid_generated_states;
            ++repeated_invalid_states;
            valid_state_generated_flag = TRUE;
            generate_new_state(user_random_generator);
            current_generated_state.cost =
                (*user_cost_function)
                (current_generated_state.parameter,
                 parameter_minimum,
                 parameter_maximum,
                 &valid_state_generated_flag,
                 OPTIONS,
                 user_data);
            if (repeated_invalid_states >
                (int) LIMIT_INVALID_GENERATED_STATES) {
                *exit_status = (int) TOO_MANY_INVALID_STATES;
                vfsr_exit(user_cost_function,
                          parameter_initial_final,
                          &final_cost,
                          exit_status);
                return (final_cost);
            }
        } while (valid_state_generated_flag == FALSE);
        --number_invalid_generated_states;
    } /* USER_INITIAL_PARAMETERS */

    /* set all states to the last one generated */
    VFOR(index_v) {
        /* ignore parameters that have too small a range */
        if (PARAMETER_RANGE_TOO_SMALL(index_v))
            continue;
        best_generated_state.parameter[index_v] =
            last_saved_state.parameter[index_v] =
            current_generated_state.parameter[index_v];
    }

    /* set all costs to the last one generated */
    best_generated_state.cost = last_saved_state.cost =
        current_generated_state.cost;

    accepted_to_generated_ratio = ONE;

    /* compute 1/number_parameters for efficiency */
    one_over_number_parameters = ONE / number_parameters;
    temperature_rescale_power = ONE / pow(
                                             (double) REANNEAL_RESCALE,
                                             one_over_number_parameters);

    tmp_var_db1 = -log(((double) TEMPERATURE_RATIO_SCALE));

    /* compute temperature scales */
    tmp_var_db2 = log((double) TEMPERATURE_ANNEAL_SCALE);
    temperature_scale = tmp_var_db1 * exp(-tmp_var_db2
                                          / (double) number_parameters);

    temperature_scale_parameters = temperature_scale;
    temperature_scale_cost = temperature_scale
        * (double) COST_PARAMETER_SCALE;

    /* do not calculate curvatures initially */
    curvature_flag = FALSE;

    /* reset the current cost and the number of generations performed */
    number_invalid_generated_states = 0;
    best_number_generated_saved =
        number_generated = recent_number_generated = 0;
    VFOR(index_v) {
        /* ignore parameters that have too small a range */
        if (PARAMETER_RANGE_TOO_SMALL(index_v))
            continue;
        index_parameter_generations[index_v] = 1;
    }

    /* MAIN ANNEALING LOOP */
    while (number_accepted <= ((LONG_INT) LIMIT_ACCEPTANCES)) {

        /* CALCULATE NEW TEMPERATURES */

        /* calculate new parameter temperatures */
        VFOR(index_v) {
            /* skip parameters with too small a range */
            if (PARAMETER_RANGE_TOO_SMALL(index_v))
                continue;

            log_new_temperature_ratio = -temperature_scale_parameters *
                pow((double) index_parameter_generations[index_v],
                    one_over_number_parameters);
            /* check (and correct) for too large an exponent */
            log_new_temperature_ratio =
                EXPONENT_CHECK(log_new_temperature_ratio);
            current_parameter_temperature[index_v] =
                initial_parameter_temperature[index_v]
                * exp(log_new_temperature_ratio);

            /* check for too small a parameter temperature */
            if (exp(log_new_temperature_ratio) < (double) SMALL_FLOAT) {
                *exit_status = (int) P_TEMP_TOO_SMALL;
                index_exit_v = index_v;

/* Note that this abridged exit can hold and print old values
   of some variables, such as current_cost_temperature. */
                vfsr_exit(user_cost_function,
                          parameter_initial_final,
                          &final_cost,
                          exit_status);
                return (final_cost);
            }
        }

        /* calculate new cost temperature */
        log_new_temperature_ratio =
            -temperature_scale_cost * pow((double)
                     index_cost_acceptances, one_over_number_parameters);
        log_new_temperature_ratio =
            EXPONENT_CHECK(log_new_temperature_ratio);
        current_cost_temperature = initial_cost_temperature
            * exp(log_new_temperature_ratio);

        /* check for too small a cost temperature */
        if (exp(log_new_temperature_ratio) < (double) SMALL_FLOAT) {
            *exit_status = (int) C_TEMP_TOO_SMALL;
/* Note that this abridged exit can hold and print old values
   of some variables. */
            vfsr_exit(user_cost_function,
                      parameter_initial_final,
                      &final_cost,
                      exit_status);
            return (final_cost);
        }
        /* GENERATE NEW PARAMETERS */

        /* generate a new valid set of parameters */
        repeated_invalid_states = 0;
        do {
            ++number_invalid_generated_states;
            ++repeated_invalid_states;
            valid_state_generated_flag = TRUE;
            generate_new_state(user_random_generator);
            current_generated_state.cost =
                (*user_cost_function)
                (current_generated_state.parameter,
                 parameter_minimum,
                 parameter_maximum,
                 &valid_state_generated_flag,
                 OPTIONS,
                 user_data);
            if (repeated_invalid_states >
                (int) LIMIT_INVALID_GENERATED_STATES) {
                *exit_status = (int) TOO_MANY_INVALID_STATES;
                vfsr_exit(user_cost_function,
                          parameter_initial_final,
                          &final_cost,
                          exit_status);
                return (final_cost);
            }
        } while (valid_state_generated_flag == FALSE);
        --number_invalid_generated_states;

        /* ACCEPT/REJECT NEW PARAMETERS */

        /* decide to accept/reject the new state */
        accept_new_state(user_random_generator);

        /* calculate the ratio of acceptances to generated states */
        accepted_to_generated_ratio =
            (double) (recent_number_acceptances + 1) /
            (double) (recent_number_generated + 1);

        /* CHECK FOR NEW MINIMUM */

        if (current_generated_state.cost < best_generated_state.cost) {
            /* NEW MINIMUM FOUND */

            fprintf (stderr, "T = %f, cost = %f, dE = %f\n",
                     current_cost_temperature, current_generated_state.cost,
                     (current_generated_state.cost - best_generated_state.cost));

            /* reset the recent acceptances and generated counts */
            recent_number_acceptances = recent_number_generated = 0;
            best_number_generated_saved = number_generated;
            best_number_accepted_saved = number_accepted;

            /* copy the current state into the best_generated state */
            best_generated_state.cost = current_generated_state.cost;
            VFOR(index_v) {
                /* ignore parameters that have too small a range */
                if (PARAMETER_RANGE_TOO_SMALL(index_v))
                    continue;
                best_generated_state.parameter[index_v] =
                    current_generated_state.parameter[index_v];
            }

            /* printout the new minimum state and value */
        }
        /* PERIODIC TESTING/REANNEALING/PRINTING SECTION */

        temp_var_int = (int) (number_accepted %
                              ((LONG_INT) TESTING_FREQUENCY_MODULUS));
        if ((temp_var_int == 0 && number_acceptances_saved
             == number_accepted)
            || accepted_to_generated_ratio
            < ((double) ACCEPTED_TO_GENERATED_RATIO)) {
            if (accepted_to_generated_ratio
                < ((double) ACCEPTED_TO_GENERATED_RATIO))
                recent_number_acceptances = recent_number_generated = 0;

            /* if best.cost repeats MAXIMUM_COST_REPEAT then exit */
            if (fabs(last_saved_state.cost
                     - best_generated_state.cost) <
                (double) COST_PRECISION) {
                ++index_cost_repeat;
                if (index_cost_repeat == ((int) MAXIMUM_COST_REPEAT)) {
                    *exit_status = (int) COST_REPEATING;
/* Note that this abridged exit can hold and print old values
   of some variables. */
                    vfsr_exit(user_cost_function,
                              parameter_initial_final,
                              &final_cost,
                              exit_status);
                    return (final_cost);
                }
            } else
                index_cost_repeat = 0;

            if (ACTIVATE_REANNEAL) {        /* set to FALSE to skip reannealing */
                /* calculate tangents, not curvatures, to reanneal */
                curvature_flag = FALSE;
                cost_derivatives(user_cost_function);
                reanneal();
            }
        }
    }

    /* FINISHED ANNEALING and MINIMIZATION */

    *exit_status = (int) NORMAL_EXIT;
    vfsr_exit(user_cost_function,
              parameter_initial_final,
              &final_cost,
              exit_status);
    return (final_cost);
}

static double vfsr_myrand (void);
static double vfsr_randflt (void);
static void vfsr_initialize_rng (void);

double vfsr_std_rng(                        /* return final_cost */
               double (*user_cost_function) (),    /* user's cost function */
               int g_number_parameters,            /* number of parameters */
               int *g_parameter_type,              /* INTEGER_TYPE, REAL_TYPE */
               double *parameter_initial_final,    /* initial and final
                                                      parameters */
               double *g_parameter_minimum,        /* minimum parameter values */
               double *g_parameter_maximum,        /* maximum parameter values */
               double *g_tangents,         /* numerical 1st derivatives of cost */
               double *g_curvature,        /* numerical 2nd derivatives of cost */
               int *exit_status,   /* exit code */
               VFSR_DEFINES * g_OPTIONS /* [optional] */,
               void *g_user_data /* passed to user_cost_function*/)
{
    vfsr_initialize_rng ();

    return vfsr (user_cost_function,
                 vfsr_randflt,
                 g_number_parameters,
                 g_parameter_type,
                 parameter_initial_final,
                 g_parameter_minimum,
                 g_parameter_maximum,
                 g_tangents,
                 g_curvature,
                 exit_status,
                 g_OPTIONS,
                 g_user_data);
}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: vfsr_exit
*        This procedures cleans up the system and frees storage, and
*        copies the best parameters and cost into final_cost and
*        parameter_initial_final
* Parameters:
* Inputs:
*        double user_cost_function - user's cost function
        int *exit_status           - exit code
* Outputs:
*        double *parameter_initial_final   - final parameters state
*        double *final_cost - final cost at state parameter_initial_final
* Global Variables:
*        curvature_flag - reset to compute curvature
*        best_generated_state.parameter    - array of best parameters
*        current_generated_state.parameter - free this storage
*        last_saved_state.parameter        - free this storage
*        best_generated_state.parameter    - free this storage
*        current_parameter_temperature     - free this storage
*        initial_parameter_temperature     - free this storage
*        index_parameter_generations       - free this storage
*        ptr_vfsr_out - output file
* Returns: None
* Calls:
*        cost_derivatives()
*        print_state()
*        free()
*        fclose()
* Description:
* Local Variables:
*        int index_v - iteration index
* Portability:
* Other:
***********************************************************************/
static void vfsr_exit(
                  double (*user_cost_function) (),        /* user's cost function */
                  double *parameter_initial_final,        /* initial and final
                                                             parameters */
                  double *final_cost,        /* best cost found to return */
                  int *exit_status           /* exit code */ )
{
    int index_v;                /* iteration index */

    /* calculate curvatures and tangents at best point */
    curvature_flag = TRUE;
    cost_derivatives(user_cost_function);

    /* return final function minimum and associated parameters */
    *final_cost = best_generated_state.cost;
    VFOR(index_v) {
        parameter_initial_final[index_v] =
            best_generated_state.parameter[index_v];
    }

    /* return unused space */
    free(current_generated_state.parameter);
    free(last_saved_state.parameter);
    free(best_generated_state.parameter);
    free(current_parameter_temperature);
    free(initial_parameter_temperature);
    free(index_parameter_generations);

    vfsr_open = FALSE;
}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: generate_new_state
* Parameters:
* Inputs:
*        double (*user_random_generator) () - the random number generator
* Outputs: None
* Global Variables:
* Returns: None
* Calls:
*        fabs                     - absolute value function
*        generate_vfsr_state      - to get a vfsr distributed state
* Description:
*        Generates a valid new state from the old state
* Local Variables:
*        int index_v               - iteration index v
*        double x                  - the new parameter value
*        double parameter_v        - the vth parameter value
*        double min_parameter_v    - the vth parameter lower bound
*        double max_parameter_v    - the vth parameter upper bound
*        double temperature_v      - the vth parameter temperature
*        double parameter_range_v  - the vth parameter range
* Portability:
* Other:
***********************************************************************/
 /* generate a new state from the old state that lies between
    the min and max parameter values */
static void generate_new_state(double (*user_random_generator) ())
{
    int index_v;
    double x;
    double parameter_v, min_parameter_v, max_parameter_v, temperature_v,
     parameter_range_v;

    /* generate a new value for each parameter */
    VFOR(index_v) {
        min_parameter_v = parameter_minimum[index_v];
        max_parameter_v = parameter_maximum[index_v];
        parameter_range_v = max_parameter_v - min_parameter_v;

        /* ignore parameters that have too small a range */
        if (fabs(parameter_range_v) < (double) SMALL_FLOAT)
            continue;

        temperature_v = current_parameter_temperature[index_v];
        parameter_v = last_saved_state.parameter[index_v];

        /* handle discrete integer parameters */
        if (INTEGER_PARAMETER(index_v)) {
            min_parameter_v -= HALF;
            max_parameter_v += HALF;
            parameter_range_v = max_parameter_v - min_parameter_v;
        }
        /* generate a new state x until within the parameter bounds */
        for (;;) {
            x = parameter_v
                + generate_vfsr_state(temperature_v, user_random_generator)
                * parameter_range_v;
            /* exit the loop if within its valid parameter range */
            if (x <= max_parameter_v - (double) SMALL_FLOAT
                && x >= min_parameter_v + (double) SMALL_FLOAT)
                break;
        }

        /* handle discrete integer parameters */
/* You might have to check rounding on your machine. */
        if (INTEGER_PARAMETER(index_v)) {
            if (x < min_parameter_v + HALF)
                x = min_parameter_v + HALF + (double) SMALL_FLOAT;
            if (x > max_parameter_v - HALF)
                x = max_parameter_v - HALF + (double) SMALL_FLOAT;

            if (x + HALF > ZERO) {
                x = (int) (x + HALF);
            } else {
                x = (int) (x - HALF);
            }
        }
        /* save the newly generated value */
        current_generated_state.parameter[index_v] = x;
    }

}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: generate_vfsr_state
* Parameters:
* Inputs: double Temperature
*        double user_random_generator() returns random number between 0,1
* Outputs: None
* Global Variables: None
* Returns: double: A random variable VFSR distributed
* Calls:
*        double pow() power function
*        double fabs() absolute value function
* Description:
*        This function generates a single value according to the
*        VFSR generating function and the passed temperature
* Local Variables:
*        double x - temporary variable
*        double y - temporary variable
*        double z - temporary variable
* Portability:
* Fully portable
* Other:
***********************************************************************/
static double generate_vfsr_state(double temp,
                           double (*user_random_generator) ())
{
    double x, y, z;

    x = (*user_random_generator) ();
    y = x < HALF ? -ONE : ONE;
    z = y * temp * (pow((ONE + ONE / temp), fabs(TWO * x - ONE)) - ONE);

    return (z);

}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: accept_new_state
*        This procedure accepts or rejects a newly generated state,
*        depending on whether the difference between new and old
*        cost functions passes a statistical test. If accepted,
*        the current state is updated.
* Parameters:
* Inputs: None
* double user_random_generator() returning a random number between 0,1
* Outputs: None
* Global Variables:
* Returns: None
* Calls:
* Description:
* Local Variables:
* Portability:
* Other:
***********************************************************************/
 /* determine whether to accept or reject the new state */
void accept_new_state(double (*user_random_generator) ())
{
    double delta_cost;
    int index_v;

    /* update accepted and generated count */
    ++number_acceptances_saved;
    ++recent_number_generated;
    ++number_generated;

    /* increment the parameter index generation for each parameter */
    VFOR(index_v) {
        /* ignore parameters with too small a range */
        if (PARAMETER_RANGE_TOO_SMALL(index_v))
            continue;
        ++(index_parameter_generations[index_v]);
    }

    /* effective cost function for testing acceptance criteria,
       calculate the cost difference and divide by the temperature */
    delta_cost = (current_generated_state.cost - last_saved_state.cost)
        / (current_cost_temperature + SMALL_FLOAT);

    /* decide to accept/reject the new state */
    if (exp(EXPONENT_CHECK(-delta_cost)) > (*user_random_generator) ()) {
        /* copy the current generated parameters to the last saved state */
        last_saved_state.cost = current_generated_state.cost;
        VFOR(index_v) {
            /* ignore parameters with too small a range */
            if (PARAMETER_RANGE_TOO_SMALL(index_v))
                continue;
            last_saved_state.parameter[index_v] =
                current_generated_state.parameter[index_v];
        }
        /* update acceptance counts */
        ++recent_number_acceptances;
        ++number_accepted;
        ++index_cost_acceptances;
        number_acceptances_saved = number_accepted;
    }
}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: reanneal
* Parameters:
* Inputs: None
* Outputs: None
* Global Variables:
* Returns: None
* Calls:
* Description:
*        Readjust temperatures of generating and acceptance functions
* Local Variables:
*        int index_v                                  - vth iteration index
*        double tmp_var_db2                           -
*        double new_temperature                       - the new temperature
*        double log_initial_current_temperature_ratio - for efficiency
* Portability:
* Other:
***********************************************************************/
 /* REANNEALING - readjust the temperatures */
static void reanneal(void)
{
    int index_v;
    double tmp_var_db2;
    double new_temperature;        /* the temperature */
    double log_initial_current_temperature_ratio;

    VFOR(index_v) {
        /* skip parameters with too small range or integer parameters */
        if (INCLUDE_INTEGER_PARAMETERS) {
            if (PARAMETER_RANGE_TOO_SMALL(index_v))
                continue;
        } else {
            if (PARAMETER_RANGE_TOO_SMALL(index_v) ||
                INTEGER_PARAMETER(index_v))
                continue;
        }

        /* ignore parameters with too small a range */
        if (PARAMETER_RANGE_TOO_SMALL(index_v))
            continue;
        /* ignore parameters with too small tangents */
        if (fabs(tangents[index_v]) < (double) SMALL_FLOAT)
            continue;

        /* reset the index of parameter generations appropriately */
        new_temperature =
            current_parameter_temperature[index_v] *
            (maximum_tangent / fabs(tangents[index_v]));
        if (new_temperature < initial_parameter_temperature[index_v]) {
            log_initial_current_temperature_ratio =
                log(((double) SMALL_FLOAT + initial_parameter_temperature[index_v])
                    / ((double) SMALL_FLOAT + new_temperature));
            index_parameter_generations[index_v] =
                (LONG_INT) ((double) SMALL_FLOAT
                            + pow(log_initial_current_temperature_ratio
                                  / temperature_scale_parameters,
                                  (double) number_parameters));
        } else
            index_parameter_generations[index_v] = 1;

        /* Reset index_parameter_generations if index reset too large,
           and also reset the initial_parameter_temperature, to achieve
           the same new temperature. */
        while (index_parameter_generations[index_v]
               > ((LONG_INT) MAXIMUM_REANNEAL_INDEX)) {
            log_new_temperature_ratio = -temperature_scale_parameters *
                pow((double) index_parameter_generations[index_v],
                    one_over_number_parameters);
            log_new_temperature_ratio =
                EXPONENT_CHECK(log_new_temperature_ratio);
            new_temperature =
                initial_parameter_temperature[index_v]
                * exp(log_new_temperature_ratio);
            index_parameter_generations[index_v] /=
                (LONG_INT) REANNEAL_RESCALE;
            initial_parameter_temperature[index_v] =
                new_temperature * pow(
                initial_parameter_temperature[index_v] / new_temperature,
                                         temperature_rescale_power);
        }
    }

    /* reanneal : Reset the index of cost acceptances to take
       into account finer detail in cost terrain. */

    /* set the starting cost_temperature to last cost found so far */
    if (initial_cost_temperature > fabs(last_saved_state.cost))
        initial_cost_temperature = fabs(last_saved_state.cost);

    if (current_cost_temperature > fabs(best_generated_state.cost)) {
        tmp_var_db2 =
            log(((double) SMALL_FLOAT + initial_cost_temperature) /
                ((double) SMALL_FLOAT + fabs(best_generated_state.cost)));
        index_cost_acceptances = (LONG_INT) ((double) SMALL_FLOAT
                                             + pow(tmp_var_db2
                                                / temperature_scale_cost,
                                            (double) number_parameters));
    } else {
        log_initial_current_temperature_ratio =
            log(((double) SMALL_FLOAT + initial_cost_temperature) /
                ((double) SMALL_FLOAT + current_cost_temperature));
        index_cost_acceptances = (LONG_INT) ((double) SMALL_FLOAT
                              + pow(log_initial_current_temperature_ratio
                                    / temperature_scale_cost,
                                    (double) number_parameters));
    }

    /* reset index_cost_temperature if index reset too large */
    while (index_cost_acceptances > ((LONG_INT) MAXIMUM_REANNEAL_INDEX)) {
        log_new_temperature_ratio = -temperature_scale_cost *
            pow((double) index_cost_acceptances,
                one_over_number_parameters);
        log_new_temperature_ratio =
            EXPONENT_CHECK(log_new_temperature_ratio);
        new_temperature = initial_cost_temperature
            * exp(log_new_temperature_ratio);
        index_cost_acceptances /= (double) REANNEAL_RESCALE;
        initial_cost_temperature =
            new_temperature * pow(
                              initial_cost_temperature / new_temperature,
                                     temperature_rescale_power);
    }
}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure: cost_derivatives
*        This procedure calculates the derivatives of the cost function
*        with respect to its parameters.  The first derivatives are
*        used as a sensitivity measure for reannealing.  The second
*        derivatives are calculated only if curvature_flag=TRUE;
*        these are a measure of the covariance of the fit when a
*        minimum is found.
* Parameters:
* Inputs: (double (*user_cost_function) () - the user's cost function
* Outputs: None
* Global Variables:
* Returns: None
* Calls:
* Description:
* Local Variables:
*        int index_v, index_vv, index_v_vv - iteration indices
*        LONG_INT saved_number_invalid_generated_states - temporary variable
*        double parameter_v, parameter_vv
*                - values of the v_th and vv_th parameters
*        double recent_best_cost - the most recently found best cost
*        double new_cost_state_1, new_cost_state_2, new_cost_state_3;
*                - cost values taken at sample points 1,2 and 3
* Portability:
* Other:
***********************************************************************/
 /* Calculate the numerical derivatives of the best
    generated state found so far */

/* In this implementation of VFSR, no checks are made for
   valid_state_generated_flag=FALSE for differential neighbors
   to the current best state. */

/* Assuming no information is given about the metric of the parameter
   space, use simple Cartesian space to calculate curvatures. */

static void cost_derivatives(double (*user_cost_function) ())
{
    int index_v, index_vv, index_v_vv, index_vv_v;
    LONG_INT saved_number_invalid_generated_states;
    double parameter_v, parameter_vv, recent_best_cost;
    double new_cost_state_1, new_cost_state_2, new_cost_state_3;

    /* save the best cost */
    recent_best_cost = best_generated_state.cost;

    /* copy the best state into the current state */
    VFOR(index_v) {
        /* ignore parameters with too small ranges */
        if (PARAMETER_RANGE_TOO_SMALL(index_v))
            continue;
        current_generated_state.parameter[index_v] =
            best_generated_state.parameter[index_v];
    }

    /* set parameters (& possibly constraints) to best state */
    saved_number_invalid_generated_states =
        number_invalid_generated_states;
    valid_state_generated_flag = TRUE;
    current_generated_state.cost =
        (*user_cost_function)
        (current_generated_state.parameter,
         parameter_minimum,
         parameter_maximum,
         &valid_state_generated_flag,
         OPTIONS,
         user_data);
    number_invalid_generated_states =
        saved_number_invalid_generated_states;
    valid_state_generated_flag = TRUE;

    VFOR(index_v) {
        /* skip parameters with too small range or integer parameters */
        if (INCLUDE_INTEGER_PARAMETERS) {
            if (PARAMETER_RANGE_TOO_SMALL(index_v)) {
                tangents[index_v] = ZERO;
                if (curvature_flag == TRUE) {
                    curvature[index_v] = ZERO;
                    continue;
                }
            }
        } else {
            if (PARAMETER_RANGE_TOO_SMALL(index_v) ||
                INTEGER_PARAMETER(index_v)) {
                tangents[index_v] = ZERO;
                if (curvature_flag == TRUE) {
                    curvature[index_v] = ZERO;
                    continue;
                }
            }
        }
        /* save parameter_v */
        parameter_v = best_generated_state.parameter[index_v];

        /* generate the first sample point and
           calculate cost at current point + DELTA_X */
        current_generated_state.parameter[index_v] =
            (ONE + (double) DELTA_X) * parameter_v;
        /* generate the first sample point */
        current_generated_state.cost =
            (*user_cost_function)
            (current_generated_state.parameter,
             parameter_minimum,
             parameter_maximum,
             &valid_state_generated_flag,
             OPTIONS,
             user_data);
        new_cost_state_1 = current_generated_state.cost;
        valid_state_generated_flag = TRUE;

        /* restore the parameter state */
        current_generated_state.parameter[index_v] = parameter_v;

        /* calculate the numerical derivative */
        tangents[index_v] = (new_cost_state_1 - recent_best_cost)
            / (((double) DELTA_X) * parameter_v + (double) SMALL_FLOAT);

        /* calculate the diagonal curvatures */
        if (curvature_flag == TRUE) {

            /* generate the second sample point and
               calculate cost at current point - DELTA_X */
            current_generated_state.parameter[index_v] =
                (ONE - (double) DELTA_X) * parameter_v;
            current_generated_state.cost =
                (*user_cost_function)
                (current_generated_state.parameter,
                 parameter_minimum,
                 parameter_maximum,
                 &valid_state_generated_flag,
                 OPTIONS,
                 user_data);
            new_cost_state_2 = current_generated_state.cost;
            valid_state_generated_flag = TRUE;

            /* restore the parameter state */
            current_generated_state.parameter[index_v] =
                parameter_v;

            /* index_v_vv: row index_v, column index_v */
            index_v_vv = ROW_COL_INDEX(index_v, index_v);

            /* calculate and store the curvature */
            curvature[index_v_vv] =
                (new_cost_state_2 - TWO * recent_best_cost
                 + new_cost_state_1) / (((double) DELTA_X)
                     * parameter_v * parameter_v + (double) SMALL_FLOAT);
        }
    }

    /* calculate off-diagonal curvatures */
    if (curvature_flag == TRUE) {
        VFOR(index_v) {
            /* save the v_th parameter */
            parameter_v = current_generated_state.parameter[index_v];

            VFOR(index_vv) {
                /* calculate only the upper diagonal */
                if (index_v <= index_vv)
                    continue;

                /* index_v_vv: row index_v, column index_vv */
                index_v_vv = ROW_COL_INDEX(index_v, index_vv);
                index_vv_v = ROW_COL_INDEX(index_vv, index_v);

                /* skip parms with too small range or integer parameters */
                if (INCLUDE_INTEGER_PARAMETERS) {
                    if (PARAMETER_RANGE_TOO_SMALL(index_v) ||
                        PARAMETER_RANGE_TOO_SMALL(index_vv)) {
                        curvature[index_vv_v] =
                            curvature[index_v_vv] = ZERO;
                        continue;
                    }
                } else {
                    if (INTEGER_PARAMETER(index_v) ||
                        INTEGER_PARAMETER(index_vv) ||
                        PARAMETER_RANGE_TOO_SMALL(index_v) ||
                        PARAMETER_RANGE_TOO_SMALL(index_vv)) {
                        curvature[index_vv_v] =
                            curvature[index_v_vv] = ZERO;
                        continue;
                    }
                }
                /* save the vv_th parameter */
                parameter_vv =
                    current_generated_state.parameter[index_vv];

                /* generate first sample point */
                current_generated_state.parameter[index_v] =
                    (ONE + (double) DELTA_X) * parameter_v;

                current_generated_state.parameter[index_vv] =
                    (ONE + (double) DELTA_X) * parameter_vv;
                current_generated_state.cost =
                    (*user_cost_function)
                    (current_generated_state.parameter,
                     parameter_minimum,
                     parameter_maximum,
                     &valid_state_generated_flag,
                     OPTIONS,
                     user_data);
                new_cost_state_1 = current_generated_state.cost;
                valid_state_generated_flag = TRUE;

                /* restore the v_th parameter */
                current_generated_state.parameter[index_v] =
                    parameter_v;

                /* generate second sample point */
                current_generated_state.cost =
                    (*user_cost_function)
                    (current_generated_state.parameter,
                     parameter_minimum,
                     parameter_maximum,
                     &valid_state_generated_flag,
                     OPTIONS,
                     user_data);
                new_cost_state_2 = current_generated_state.cost;
                valid_state_generated_flag = TRUE;

                /* restore the vv_th parameter */
                current_generated_state.parameter[index_vv] =
                    parameter_vv;

                /* generate third sample point */
                current_generated_state.parameter[index_v] =
                    (ONE + (double) DELTA_X)
                    * best_generated_state.parameter[index_v];
                current_generated_state.cost =
                    (*user_cost_function)
                    (current_generated_state.parameter,
                     parameter_minimum,
                     parameter_maximum,
                     &valid_state_generated_flag,
                     OPTIONS,
                     user_data);
                new_cost_state_3 = current_generated_state.cost;
                valid_state_generated_flag = TRUE;

                /* restore the v_th parameter */
                current_generated_state.parameter[index_v] =
                    parameter_v;

                /* calculate and store the curvature */
                curvature[index_vv_v] = curvature[index_v_vv] =
                    (new_cost_state_1 - new_cost_state_2
                     - new_cost_state_3 + recent_best_cost)
                    / (((double) DELTA_X) * ((double) DELTA_X)
                       * parameter_v * parameter_vv
                       + (double) SMALL_FLOAT);
            }
        }
    }
    /* restore the best cost function value */
    current_generated_state.cost = recent_best_cost;

    /* find the maximum |tangent| from all tangents */
    maximum_tangent = 0;
    VFOR(index_v) {
        /* ignore too small ranges and integer parameters types */
        if (INCLUDE_INTEGER_PARAMETERS) {
            if (PARAMETER_RANGE_TOO_SMALL(index_v))
                continue;
        } else {
            if (PARAMETER_RANGE_TOO_SMALL(index_v) ||
                INTEGER_PARAMETER(index_v))
                continue;
        }

        /* find the maximum |tangent| (from all tangents) */
        if (fabs(tangents[index_v]) > maximum_tangent) {
            maximum_tangent = fabs(tangents[index_v]);
        }
    }

}

#define SHUFFLE 256                /* size of random array */
#define MULT 25173
#define MOD ((long int) 65536)
#define INCR 13849
#define FMOD ((double) 65536.0)

static long int seed = 696969;
double random_array[SHUFFLE];        /* random variables */

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure:
*        double myrand(void) - returns a random number between 0 and 1
* Parameters: None
* Inputs:
* Outputs:
* Global Variables:
*        double seed             - The rng seed
* Returns: None
* Calls: None
* Description:
*        This routine returns the random number generator between 0 and 1
* Local Variables: None
* Portability:
* Other:
***********************************************************************/

static double vfsr_myrand(void)
/* returns random number in {0,1} */
{
    seed = (MULT * seed + INCR) % MOD;
    return ((double) seed / FMOD);
}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure:
*        double randflt(void) - Shuffles random numbers in random_array[]
* Parameters: None
* Inputs:
* Outputs:
* Global Variables: None
* Returns: None
* Calls:
*        myrand()
* Description:
*        This routine initializes the random number generator
* Local Variables: None
* Portability:
* Other:
***********************************************************************/

static double vfsr_randflt(void)
/* shuffles random numbers in random_array[SHUFFLE] array */
{
    double rranf;
    int kranf;

    kranf = (int) (vfsr_myrand() * SHUFFLE) % SHUFFLE;
    rranf = *(random_array + kranf);
    *(random_array + kranf) = vfsr_myrand();

    return (rranf);
}

/***********************************************************************
* Author: Lester Ingber, Bruce Rosen (copyright) (c)
* Date         5 Nov 92
* Procedure:
*        initialize_rng() - to initialize the random number generator
* Parameters: None
* Inputs:
* Outputs:
* Global Variables: None
* Returns: None
* Calls:
*        myrand()
*        randflt()
* Description:
*        This routine initializes the random number generator
* Local Variables: None
* Portability:
* Other:
***********************************************************************/

static void vfsr_initialize_rng(void)
{
    int n;
    double x;

    for (n = 0; n < SHUFFLE; ++n)
        random_array[n] = vfsr_myrand();
    for (n = 0; n < 1000; ++n)        /* warm up random generator */
        x = vfsr_randflt();
}
