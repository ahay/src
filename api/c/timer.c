/* Execution timer */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include <sys/time.h>
#include "timer.h"
#include "alloc.h"

#include "_defs.h"
#include "_bool.h"
/*^*/

#ifndef _sf_timer_h

typedef struct ExecTimer *sf_timer;
/* abstract data type */
/*^*/

#endif

struct ExecTimer {
    struct timeval start_time;
    /* Time difference between the last start and stop */
    float  diff_time;
    /* TOTAL time difference between starts and stops */
    float  total_time;
    /* flag if the stop watch is running */
    bool running;
    /* Number of times clock has been started
       and stopped to allow averaging */
    int clock_sessions;
};
/* concrete data type */

static float sf_get_diff_time (sf_timer timer) {
    struct timeval t_time;
    gettimeofday (&t_time, 0);

    /* time difference in msec. */
    return (float)(1000.0*(t_time.tv_sec - timer->start_time.tv_sec) + 
                   (0.001*(t_time.tv_usec - timer->start_time.tv_usec)));
}

sf_timer sf_timer_init (void)
/*< Initialize timer object. >*/
{
    sf_timer timer;

    timer = (sf_timer)sf_alloc (1, sizeof (struct ExecTimer));

    timer->diff_time = 0.0;
    timer->total_time = 0.0;
    timer->running = false;
    timer->clock_sessions = 0;
    return timer;
}

void sf_timer_close (sf_timer timer)
/*< Destroy timer object. >*/
{
    free (timer);
}

void sf_timer_start (sf_timer timer)
/*< Start time measurement session. >*/
{
    gettimeofday (&timer->start_time, 0);
    timer->running = true;
}

void sf_timer_stop (sf_timer timer)
/*< Stop time measurement session. >*/
{
    timer->diff_time = sf_get_diff_time (timer);
    timer->total_time += timer->diff_time;
    timer->running = false;
    timer->clock_sessions++;
}

void sf_timer_reset (sf_timer timer)
/*< Reset the timer to 0. Does not change the timer running state. >*/
{
    timer->diff_time = 0;
    timer->total_time = 0;
    timer->clock_sessions = 0;
    if (timer->running)
        gettimeofday (&timer->start_time, 0);
}

float sf_timer_get_total_time (sf_timer timer) 
/*< Total time in msec for all sessions. after start. >*/
{
    float retval = timer->total_time;
    if (timer->running)
        retval += sf_get_diff_time (timer);

    return retval;
}

float sf_timer_get_diff_time (sf_timer timer) 
/*< Time in msec for the last session. >*/
{
    return (timer->running) ? sf_get_diff_time (timer)
                            : timer->diff_time;
}

float sf_timer_get_average_time (sf_timer timer)
/*< Average time in msec. for all runs. >*/
{
    return timer->total_time/timer->clock_sessions;
}

