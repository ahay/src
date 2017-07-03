/* Point in 5-D phase-space grid corresponding to 3-D medium */
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

#ifndef _esc_point3_h

#ifdef sun
#include <inttypes.h>
#else
#include <stdint.h>
#endif
/*^*/

typedef struct EscPoint3 *sf_esc_point3;
/* abstract data type */
/*^*/

/* Colors for different boundaries */
typedef enum { ESC3_LEFT = 1 << 0, ESC3_RIGHT = 1 << 1,
               ESC3_TOP = 1 << 2, ESC3_BOTTOM = 1 << 3,
               ESC3_NEAR = 1 << 4, ESC3_FAR = 1 << 5,
               ESC3_INSIDE = 1 << 6 } EscColor3;
typedef enum { ESC3_AXIS_Z = 0, ESC3_AXIS_X = 1, ESC3_AXIS_Y = 2,
               ESC3_AXIS_A = 3, ESC3_AXIS_B = 4, ESC3_DIMS = 5 } EscAxisDim3;
/* Direction along an axis */
typedef enum { ESC3_BACK = 0, ESC3_FORW = 1, ESC3_DIRS = 2 } EscDirection3;
/* Escape variables */
#ifdef ESC_EQ_WITH_L
typedef enum { ESC3_Z = 0, ESC3_X = 1, ESC_Y = 2, ESC3_T = 3, ESC3_L = 4,
               ESC3_NUM = 5 } EscType3;
#else
typedef enum { ESC3_Z = 0, ESC3_X = 1, ESC3_Y = 2, ESC3_T = 3,
               ESC3_NUM = 4 } EscType3;
#endif
/*^*/

extern const char* sf_esc_point3_str[ESC3_NUM];
/*^*/

#endif

#ifdef ESC_EQ_WITH_L
const char* sf_esc_point3_str[ESC3_NUM] = { "z", "x", "y", "t", "l" };
#else
const char* sf_esc_point3_str[ESC3_NUM] = { "z", "x", "y", "t" };
#endif

struct EscPoint3 {
    float    e[ESC3_NUM]; /* Escape variables */
    uint32_t type:1; /* Child or parent */
    uint32_t nc:4; /* Number of connections */
    uint32_t f:10; /* Connections in all directions */
    uint32_t col:6; /* Color (origin) */
};
/* concrete data type */

int sf_esc_point3_sizeof (void)
/*< Returns size of object in bytes >*/
{
    return sizeof (struct EscPoint3);
}

void sf_esc_point3_reset (sf_esc_point3 esc_point)
/*< Reset object to default state >*/
{
    int i;

    for (i = 0; i < ESC3_NUM; i++)
        esc_point->e[i] = 0.0;

    esc_point->type = 0;
    esc_point->nc = 0;
    esc_point->f = 0;
    esc_point->col = 0;
}

sf_esc_point3 sf_esc_point3_init (void)
/*< Initialize object >*/
{
    sf_esc_point3 esc_point = (sf_esc_point3)sf_alloc (1, sizeof (struct EscPoint3));

    sf_esc_point3_reset (esc_point);

    return esc_point;
}

void sf_esc_point3_close (sf_esc_point3 esc_point)
/*< Destroy object >*/
{
    free (esc_point);
}

float sf_esc_point3_get_esc_var (sf_esc_point3 esc_point, EscType3 i)
/*< Get escape length >*/
{
    return esc_point->e[i];
}

EscColor3 sf_esc_point3_get_col (sf_esc_point3 esc_point)
/*< Get point color >*/
{
    return esc_point->col;
}

void sf_esc_point3_set_esc_var (sf_esc_point3 esc_point, EscType3 i, float f)
/*< Set escape variable >*/
{
    esc_point->e[i] = f;
}

void sf_esc_point3_set_col (sf_esc_point3 esc_point, EscColor3 col)
/*< Set point color >*/
{
    esc_point->col = col;
}

bool sf_esc_point3_is_child (sf_esc_point3 esc_point)
/*< Return true, if point is a parent >*/
{
    return (0 == esc_point->type);
}

bool sf_esc_point3_is_parent (sf_esc_point3 esc_point)
/*< Return true, if point is a parent >*/
{
    return (1 == esc_point->type);
}

void sf_esc_point3_become_parent (sf_esc_point3 esc_point)
/*< Change type to parent >*/
{
    if (1 == esc_point->type)
        sf_error ("sf_esc_point3: point is already a parent");

    esc_point->type = 1;
}

void sf_esc_point3_add_parent_link (sf_esc_point3 esc_point, EscAxisDim3 dim,
                                    EscDirection3 dir)
/*< Remove parent in direction (dir) along axis (dim) >*/
{
    if (esc_point->f & (1 << (2*dim)))
        sf_error ("sf_esc_point3: parent link already exists");

    esc_point->f |= (1 << (2*dim));
    esc_point->nc++;

    if (ESC3_FORW == dir)
        esc_point->f |= (1 << (2*dim + 1));
}

void sf_esc_point3_remove_parent_link (sf_esc_point3 esc_point, EscAxisDim3 dim)
/*< Remove parent in direction (dir) >*/
{
    if (!(esc_point->f & (1 << (2*dim))))
        sf_error ("sf_esc_point3: removing non-existent parent link");

    esc_point->f ^= (1 << (2*dim));
    esc_point->nc--;

    if (esc_point->f &= (1 << (2*dim + 1)))
        esc_point->f ^= (1 << (2*dim + 1));
}

bool sf_esc_point3_has_parent_link (sf_esc_point3 esc_point, EscAxisDim3 dim,
                                    EscDirection3 *dir)
/*< Return true, if the point has a parent along axis (dim);
    dir will be set to its direction >*/
{
    if (!(esc_point->f & (1 << (2*dim))))
        return false;
    *dir = (esc_point->f & (1 << (2*dim + 1))) ? ESC3_FORW
                                               : ESC3_BACK;
    return true;
}

