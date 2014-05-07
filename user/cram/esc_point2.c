/* Point in 3-D phase-space grid corresponding to 2-D medium */
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

#ifndef _esc_point2_h

#ifdef sun
#include <inttypes.h>
#else
#include <stdint.h>
#endif
/*^*/

typedef struct EscPoint2 *sf_esc_point2;
/* abstract data type */
/*^*/

/* Colors for different boundaries */
typedef enum { ESC2_LEFT = 1 << 0, ESC2_RIGHT = 1 << 1,
               ESC2_TOP = 1 << 2, ESC2_BOTTOM = 1 << 3 } EscColor2;
typedef enum { ESC2_AXIS_Z = 0, ESC2_AXIS_X = 1, ESC2_AXIS_A = 2,
               ESC2_DIMS = 3 } EscAxisDim2;
/* Direction along an axis */
typedef enum { ESC2_BACK = 0, ESC2_FORW = 1, ESC2_DIRS = 2 } EscDirection2;
/* Escape variables */
#ifdef ESC_EQ_WITH_L
typedef enum { ESC2_Z = 0, ESC2_X = 1, ESC2_T = 2, ESC2_L = 3,
               ESC2_NUM = 4 } EscType2;
#else
typedef enum { ESC2_Z = 0, ESC2_X = 1, ESC2_T = 2,
               ESC2_NUM = 3 } EscType2;
#endif
/*^*/

extern const char* sf_esc_point2_str[ESC2_NUM];
/*^*/

#endif

const char* sf_esc_point2_str[ESC2_NUM] = { "z", "x", "t" };

struct EscPoint2 {
    float    e[ESC2_NUM]; /* Escape variables */
    uint32_t type:1; /* Child or parent */
    uint32_t nc:3; /* Number of connections */
    uint32_t f:6; /* Connections in all directions */
    uint32_t col:4; /* Color (origin) */
    uint32_t trc:4; /* Is traced */
};
/* concrete data type */

int sf_esc_point2_sizeof (void)
/*< Returns size of object in bytes >*/
{
    return sizeof (struct EscPoint2);
}

void sf_esc_point2_reset (sf_esc_point2 esc_point)
/*< Reset object to default state >*/
{
    int i;

    for (i = 0; i < ESC2_NUM; i++)
        esc_point->e[i] = 0.0;

    esc_point->type = 0;
    esc_point->nc = 0;
    esc_point->f = 0;
    esc_point->col = 0;
    esc_point->trc = 0;
}

sf_esc_point2 sf_esc_point2_init (void)
/*< Initialize object >*/
{
    sf_esc_point2 esc_point = (sf_esc_point2)sf_alloc (1, sizeof (struct EscPoint2));

    sf_esc_point2_reset (esc_point);

    return esc_point;
}

void sf_esc_point2_close (sf_esc_point2 esc_point)
/*< Destroy object >*/
{
    free (esc_point);
}

float sf_esc_point2_get_esc_var (sf_esc_point2 esc_point, EscType2 i)
/*< Get escape length >*/
{
    return esc_point->e[i];
}

EscColor2 sf_esc_point2_get_col (sf_esc_point2 esc_point)
/*< Get point color >*/
{
    return esc_point->col;
}

void sf_esc_point2_set_esc_var (sf_esc_point2 esc_point, EscType2 i, float f)
/*< Set escape variable >*/
{
    esc_point->e[i] = f;
}

void sf_esc_point2_set_col (sf_esc_point2 esc_point, EscColor2 col)
/*< Set point color >*/
{
    esc_point->col = col;
}

bool sf_esc_point2_is_child (sf_esc_point2 esc_point)
/*< Return true, if point is a parent >*/
{
    return (0 == esc_point->type);
}

bool sf_esc_point2_is_parent (sf_esc_point2 esc_point)
/*< Return true, if point is a parent >*/
{
    return (1 == esc_point->type);
}

void sf_esc_point2_become_parent (sf_esc_point2 esc_point)
/*< Change type to parent >*/
{
    if (1 == esc_point->type)
        sf_error ("sf_esc_point2: point is already a parent");

    esc_point->type = 1;
}

void sf_esc_point2_add_parent_link (sf_esc_point2 esc_point, EscAxisDim2 dim,
                                    EscDirection2 dir)
/*< Remove parent in direction (dir) along axis (dim) >*/
{
    if (esc_point->f & (1 << (2*dim)))
        sf_error ("sf_esc_point2: parent link already exists");

    esc_point->f |= (1 << (2*dim));
    esc_point->nc++;

    if (ESC2_FORW == dir)
        esc_point->f |= (1 << (2*dim + 1));
}

void sf_esc_point2_remove_parent_link (sf_esc_point2 esc_point, EscAxisDim2 dim)
/*< Remove parent in direction (dir) >*/
{
    if (!(esc_point->f & (1 << (2*dim))))
        sf_error ("sf_esc_point2: removing non-existent parent link");

    esc_point->f ^= (1 << (2*dim));
    esc_point->nc--;

    if (esc_point->f &= (1 << (2*dim + 1)))
        esc_point->f ^= (1 << (2*dim + 1));
}

bool sf_esc_point2_has_parent_link (sf_esc_point2 esc_point, EscAxisDim2 dim,
                                    EscDirection2 *dir)
/*< Return true, if the point has a parent along axis (dim);
    dir will be set to its direction >*/
{
    if (!(esc_point->f & (1 << (2*dim))))
        return false;
    *dir = (esc_point->f & (1 << (2*dim + 1))) ? ESC2_FORW
                                               : ESC2_BACK;
    return true;
}

void sf_esc_point2_set_traced (sf_esc_point2 esc_point, bool traced)
/*< Set ray traced flag >*/
{
    esc_point->trc = traced;
}

bool sf_esc_point2_is_traced (sf_esc_point2 esc_point)
/*< Return true, if point is ray traced >*/
{
    return (1 == esc_point->trc);
}

