/* Survey information (source and receiver locations) for 2-D angle migration */
/*
  Copyright (C) 2011 University of Texas at Austin

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

#include <rsf.h>

#ifndef _cram_survey2_h

typedef struct CRAMSurvey2 *sf_cram_survey2;
/* abstract data type */
/*^*/

#endif

struct CRAMSurvey2 {
    int   ns, nh;
    float s0, h0;
    float ds, dh;
    bool  absoff;
};
/* concrete data type */

sf_cram_survey2 sf_cram_survey2_init (sf_file data)
/*< Initialize object >*/
{
    bool kmah;
    sf_cram_survey2 cram_survey = (sf_cram_survey2)sf_alloc (1, sizeof (struct CRAMSurvey2));

    /* Absolute or relative (with respect to source) receiver positions */
    if (!sf_histbool (data, "absoff", &cram_survey->absoff)) sf_error ("No absoff= in data");
    /* Use KMAH phase shifts */
    if (!sf_histbool (data, "KMAH", &kmah)) sf_error ("No KMAH= in data");

    if (cram_survey->absoff)
        sf_warning ("Assuming absolute offset in data");

    if (kmah) {
        if (!sf_histint (data, "n3", &cram_survey->nh)) sf_error ("No n3= in data");
        if (!sf_histfloat (data, "o3", &cram_survey->h0)) sf_error ("No o2= in data");
        if (!sf_histfloat (data, "d3", &cram_survey->dh)) sf_error ("No d2= in data");
        if (!sf_histint (data, "n4", &cram_survey->ns)) sf_error ("No n4= in data");
        if (!sf_histfloat (data, "o4", &cram_survey->s0)) sf_error ("No o3= in data");
        if (!sf_histfloat (data, "d4", &cram_survey->ds)) sf_error ("No d3= in data");
    } else {
        if (!sf_histint (data, "n2", &cram_survey->nh)) sf_error ("No n2= in data");
        if (!sf_histfloat (data, "o2", &cram_survey->h0)) sf_error ("No o2= in data");
        if (!sf_histfloat (data, "d2", &cram_survey->dh)) sf_error ("No d2= in data");
        if (!sf_histint (data, "n3", &cram_survey->ns)) sf_error ("No n3= in data");
        if (!sf_histfloat (data, "o3", &cram_survey->s0)) sf_error ("No o3= in data");
        if (!sf_histfloat (data, "d3", &cram_survey->ds)) sf_error ("No d3= in data");
    }

    return cram_survey;
}

void sf_cram_survey2_close (sf_cram_survey2 cram_survey)
/*< Destroy object >*/
{
    free (cram_survey);
}

float sf_cram_survey2_get_src_sampling (sf_cram_survey2 cram_survey)
/*< Return source spatial sampling >*/
{
    return cram_survey->ds;
}

float sf_cram_survey2_get_rcv_sampling (sf_cram_survey2 cram_survey)
/*< Return receiver spatial sampling >*/
{
    return cram_survey->dh;
}

float sf_cram_survey2_get_first_source (sf_cram_survey2 cram_survey, int *is, int *ns)
/*< Initialize iterator through sources, return source index, number of source >*/
{
    *is = 0;
    *ns = cram_survey->ns;

    return cram_survey->s0;
}

float sf_cram_survey2_get_next_source (sf_cram_survey2 cram_survey, int *is)
/*< Iterator through sources, return source index >*/
{
    *is = *is + 1;
    if (*is >= cram_survey->ns) {
        *is = -1;
        return SF_HUGE;
    }

    return cram_survey->s0 + (*is)*cram_survey->ds;
}

float sf_cram_survey2_get_first_receiver (sf_cram_survey2 cram_survey, int is,
                                          int *ih, int *nh, size_t *i)
/*< Initialize iterator through receivers for a given source,
    return receiver index, number of receivers, and trace index >*/
{
    *ih = 0;
    *i = is*cram_survey->nh;
    *nh = cram_survey->nh;

    return (cram_survey->absoff)
           ? cram_survey->h0
           : cram_survey->s0 + is*cram_survey->ds +
             cram_survey->h0;
}

float sf_cram_survey2_get_next_receiver (sf_cram_survey2 cram_survey, int is,
                                         int *ih, size_t *i)
/*< Iterator through receivers for a given source,
    return receiver index and trace index >*/
{
    *ih = *ih + 1;
    if (*ih >= cram_survey->nh) {
        *ih = -1;
        return SF_HUGE;
    }
    *i = is*cram_survey->nh + (*ih);

    return (cram_survey->absoff)
           ? cram_survey->h0 + (*ih)*cram_survey->dh
           : cram_survey->s0 + is*cram_survey->ds +
             cram_survey->h0 + (*ih)*cram_survey->dh;
}

