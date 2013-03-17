/* Prepare survey info for 3-D angle-domain migration. */
/*
  Copyright (C) 2012 University of Texas at Austin

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

#include "cram_helper.h"

/* Comparison routine for qsort() by x values */
static int sicomp (const void *v1, const void *v2) {
    float f1 = *(float*)v1, f2 = *(float*)v2;
    if (f1 < f2)
        return -1;
    else if (f1 > f2)
        return 1;
    else
        return 0;
}

/* Comparison routine for qsort() by indices */
static int vcomp (const void *v1, const void *v2) {
    int *p1 = (int*)v1, *p2 = (int*)v2;
    int i1 = SF_MIN ((*p1), SF_MIN ((*(p1 + 1)), (*(p1 + 2))));
    int i2 = SF_MIN ((*p2), SF_MIN ((*(p2 + 1)), (*(p2 + 2))));
    if (i1 < i2)
        return -1;
    else if (i1 > i2)
        return 1;
    else
        return 0;
}

/* Do 2-D triangulation of an array of sources/recievers, return indices to
   vertices in all triangles */
int* sf_prepare_cram_survey3_tri (int n, float *s, float emax, int *ntr) {
    int i, j, *v;
    float *si, x0, x1, x2, y0, y1, y2;
    float l1, l2, l3;
    union {
        float f;
        int i;
    } uf2i;

    *ntr = 0;

    if (n <= 2)
        return NULL;

    si = sf_floatalloc (3*(n + 3));
    v = sf_intalloc (3*3*n);

    /* Copy coordinates into a separate array for sorting,
       remember their original indices */
    for (i = 0; i < n; i++) {
        si[i*3] = s[i*2];
        si[i*3 + 1] = s[i*2 + 1];
        uf2i.i = i;
        si[i*3 + 2] = uf2i.f;
    }
    /* Sort by x values for the triangulation routine */
    qsort (si, n, 3*sizeof(float), sicomp);
    /* Triangulation itself */
    if (0 != sf_cram_triangulate (n, 3, si, v, ntr))
        sf_error ("Triangulation error");
    /* Convert triangle vertex indices to trace indices */
    j = 0;
    for (i = 0; i < (*ntr); i++) {
        /* Vertices */
        x0 = si[v[i*3]*3]; y0 = si[v[i*3]*3 + 1];
        x1 = si[v[i*3 + 1]*3]; y1 = si[v[i*3 + 1]*3 + 1];
        x2 = si[v[i*3 + 2]*3]; y2 = si[v[i*3 + 2]*3 + 1];
        /* Edge length */
        l1 = hypotf (x1 - x0, y1 - y0);
        l2 = hypotf (x2 - x0, y2 - y0);
        l3 = hypotf (x2 - x1, y2 - y1);
        if (l1 <= emax && l2 <= emax && l3 <= emax) {
            /* Translate triangulation indices into trace indices */
            uf2i.f = si[v[i*3]*3 + 2];
            v[j*3] = uf2i.i;
            uf2i.f = si[v[i*3 + 1]*3 + 2];
            v[j*3 + 1] = uf2i.i;
            uf2i.f = si[v[i*3 + 2]*3 + 2];
            v[j*3 + 2] = uf2i.i;
            j++;
        }
    }
    *ntr = j;
    /* Sort by trace indices */
    qsort (v, *ntr, 3*sizeof(int), vcomp);

    free (si);

    return v;
}

int main (int argc, char* argv[]) {
    int nt, nst, *sh = NULL, **hh = NULL;
    size_t i, is, n, ns, nh, nht = 0, *irec = NULL, *ihtr = NULL, sz;
    float gx, gy, esmax, ehmax, *buf, *src, *rec, *minmax;
    sf_complex s, ps, r;
    sf_file data, survey, sxsy = NULL, gxgy = NULL;

    bool verb, tri;

    sf_init (argc, argv);

    data = sf_input ("in");
    /* Common-shot 2-D data */
    survey = sf_output ("out");
    /* Survey info */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool ("tri", &tri)) tri = false;
    /* triangulation flag */

    /* Data dimensions */
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in data");
    n = sf_leftsize (data, 1);

    if (sf_getstring ("sxsy")) {
    /* File with shot coordinates */
        sxsy = sf_input ("sxsy");
        if (SF_COMPLEX != sf_gettype (sxsy)) sf_error ("Need complex type in sxsy=");
        if (!sf_histlargeint (sxsy, "n2", (off_t*)&i)) sf_error ("No n2= in sxsy");
        if (i != n) sf_error ("n2 in sxsy= should be equal to the number of input traces");
    } else
        sf_error ("Need sxsy=");

    if (sf_getstring ("gxgy")) {
    /* File with receiver coordinates */
        gxgy = sf_input ("gxgy");
        if (SF_COMPLEX != sf_gettype (gxgy)) sf_error ("Need complex type in gxgy=");
        if (!sf_histlargeint (gxgy, "n2", (off_t*)&i)) sf_error ("No n2= in gxgy");
        if (i != n) sf_error ("n2 in gxgy= should be equal to the number of input traces");
    } else
        sf_error ("Need gxgy=");

    if (!sf_getfloat ("esmax", &esmax)) esmax = 0.2;
    /* Maximum edge length in the shot triangulation */
    if (!sf_getfloat ("ehmax", &ehmax)) ehmax = 0.1;
    /* Maximum edge length in the receiver triangulation */

    buf = sf_floatalloc (n*6);

    src = sf_floatalloc ((size_t)n*(size_t)2);
    irec = (size_t*)sf_alloc ((size_t)n, sizeof(size_t));
    rec = sf_floatalloc ((size_t)n*(size_t)2);
    minmax = sf_floatalloc ((size_t)n*(size_t)4);

    ps = sf_cmplx (SF_HUGE, SF_HUGE);
    ns = 0;
    nh = 0;
    /* Loop over source and receiver coordinates,
       assume data sorted by shots */
    while (nh < n) {
        /* Source coordinate */
        sf_complexread (&s, 1, sxsy);
        /* Receiver coordinate */
        sf_complexread (&r, 1, gxgy);
        gx = crealf (r);
        gy = cimagf (r);
        if (crealf (s) != crealf (ps) ||
            cimagf (s) != cimagf (ps)) {
            /* Next source */
            irec[ns] = nh;
            src[ns*(size_t)2] = crealf (s);
            src[ns*(size_t)2 + (size_t)1] = cimagf (s);
            ns++;
            ps = s;
            minmax[(ns - (size_t)1)*(size_t)4] = gx;
            minmax[(ns - (size_t)1)*(size_t)4 + (size_t)1] = gx;
            minmax[(ns - (size_t)1)*(size_t)4 + (size_t)2] = gy;
            minmax[(ns - (size_t)1)*(size_t)4 + (size_t)3] = gy;
        } else {
            if (gx < minmax[(ns - (size_t)1)*(size_t)4] || 0 == nh)
                minmax[(ns - (size_t)1)*(size_t)4] = gx;
            if (gx > minmax[(ns - (size_t)1)*(size_t)4 + (size_t)1] || 0 == nh)
                minmax[(ns - (size_t)1)*(size_t)4 + (size_t)1] = gx;
            if (gy < minmax[(ns - (size_t)1)*(size_t)4 + (size_t)2] || 0 == nh)
                minmax[(ns - (size_t)1)*(size_t)4 + (size_t)2] = gy;
            if (gy > minmax[(ns - (size_t)1)*(size_t)4 + (size_t)3] || 0 == nh)
                minmax[(ns - (size_t)1)*(size_t)4 + (size_t)3] = gy;
        }
        /* Next receiver */
        rec[nh*(size_t)2] = gx;
        rec[nh*(size_t)2 + (size_t)1] = gy;
        nh++;
    }

    if (verb)
        sf_warning ("%d shots found, %d receivers total", ns, nh);

    /* Total size in bytes for shot coordinates: number of receivers in each shot +
       shot coordinates + min/max coordinates for every shot + receiver coordinates */
    sz = (size_t)ns*(size_t)sizeof(size_t) +
         (size_t)ns*(size_t)(sizeof(float)*2) +
         (size_t)ns*(size_t)(sizeof(float)*4) +
         (size_t)nh*(size_t)(sizeof(float)*2);

    if (tri) {
        sh = sf_prepare_cram_survey3_tri (ns, src, esmax, &nst);
        if (verb)
            sf_warning ("%d triangles in the source triangulation", nst);
        sz += (size_t)nst*(size_t)(3*sizeof(int));

        /* Loop over shots, triangulate receivers */
        hh = (int**)sf_alloc (ns, sizeof(int*));
        ihtr = (size_t*)sf_alloc (ns, sizeof(size_t));
        sz += (size_t)ns*(size_t)sizeof(size_t);
        nht = 0;
        if (verb)
            sf_warning ("Triangulating receivers");
        for (is = 0; is < ns; is++) {
            if (verb)
                sf_warning ("Triangulating receivers in shot %d of %d;", is + 1, ns);
            hh[is] = sf_prepare_cram_survey3_tri (is != (ns - 1) ? irec[is + 1] - irec[is]
                                                                 : nh - irec[is],
                                                  &rec[irec[is]*2], ehmax, &nt);
            ihtr[is] = nht;
            nht += (size_t)nt;
            /* Size of this receiver triangulation + number of triangles in bytes */
            sz += (size_t)nt*(size_t)(3*sizeof(int));
        }
        if (verb) {
            sf_warning (".");
            sf_warning ("%lu triangles in the receiver triangulation", nht);
        }
    }

    sf_settype (survey, SF_UCHAR);
    sf_putlargeint (survey, "n1", sz);
    sf_putfloat (survey, "o1", 0.0);
    sf_putfloat (survey, "d1", 1.0);
    sf_putstring (survey, "label1", "");
    sf_putstring (survey, "unit1", "");
    sf_putint (survey, "n2", 1);
    sf_putint (survey, "n3", 1);
    sf_putint (survey, "n4", 1);
    sf_putint (survey, "n5", 1);
    sf_putint (survey, "n6", 1);

    sf_putlargeint (survey, "Nshot", ns);
    sf_putlargeint (survey, "Nrec", nh);
    sf_putstring (survey, "Istri", tri ? "y" : "n");
    if (tri) {
        sf_putint (survey, "Nstri", nst);
        sf_putlargeint (survey, "Nhtri", nht);
    }

    /* Array of number of receivers in each shot */
    sf_ucharwrite ((unsigned char*)irec, (size_t)ns*(size_t)sizeof(size_t), survey);
    /* Shot coordinates */
    sf_ucharwrite ((unsigned char*)src, (size_t)ns*(size_t)(sizeof(float)*2), survey);
    /* Min/max coordinates for every shot */
    sf_ucharwrite ((unsigned char*)minmax, (size_t)ns*(size_t)(sizeof(float)*4), survey);
    /* Receiver coordinates */
    sf_ucharwrite ((unsigned char*)rec, (size_t)nh*(size_t)(sizeof(float)*2), survey);
    if (tri) {
        /* Number of triangles in each receiver triangulation */
        sf_ucharwrite ((unsigned char*)ihtr, (size_t)ns*(size_t)sizeof(size_t), survey);
        /* Shot triangulation */
        if (sh) {
            sf_ucharwrite ((unsigned char*)sh, (size_t)nst*(size_t)(3*sizeof(int)), survey);
            free (sh);
        }
        for (is = 0; is < ns; is++) {
            /* Receiver triangulation */
            if (hh[is]) {
                if (is != (ns - 1))
                    sz = (size_t)ihtr[is + 1] - (size_t)ihtr[is];
                else
                    sz = nht - (size_t)ihtr[is];
                sf_ucharwrite ((unsigned char*)hh[is], nt*(size_t)(3*sizeof(int)), survey);
                free (hh[is]);
            }
        }
        free (hh);
        free (ihtr);
    }

    sf_fileclose (sxsy);
    sf_fileclose (gxgy);
    free (buf);
    free (irec);
    free (src);
    free (rec);

    return 0;
}

