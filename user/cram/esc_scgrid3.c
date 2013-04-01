/* 5-D phase-space escape values grid consisting of supercells */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <sys/select.h>
#include <unistd.h>
#include <netinet/in.h>
#include <netinet/tcp.h>

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

#include <rsf.h>

#ifndef _esc_scgrid3_h

#include "esc_tracer3.h"
/*^*/

typedef struct EscSCgrid3 *sf_esc_scgrid3;
/* abstract data type */
/*^*/

typedef struct {
    size_t id;          /* Message ID */
    int    ud1, ud2;    /* User defined IDs */
    int    iab;         /* Requested angle index */
    float  z, x, y;     /* Position in space */
} sf_esc_scgrid3_areq;
/* Structure for requesting one (z,x,y) point in angle space */
/*^*/

typedef struct {
    size_t id;                 /* Message ID */
    int    ud1, ud2;           /* User defined IDs */
    float  vals[ESC3_NUM + 4]; /* Position in space */
} sf_esc_scgrid3_avals;
/* Structure for getting requested escape values back */
/*^*/

#define SCGRID3_MAX_STENCIL 36
/*^*/

#endif

typedef struct {
    int    iab;            /* Initial angle index */
    float  a, b, fa, fb;   /* Azimuth, inclination, fractions */
    float  vals[ESC3_NUM]; /* Position in space */
} sf_esc_scgrid3_invals;
/* Structure for storing values in progress */

#include "einspline.h"
#include "esc_helper.h"
/*
#ifdef HAVE_MKL
#include <mkl.h>
#define MKL_LU
typedef MKL_INT lu_int;
#else
*/
typedef int lu_int;
/*
#endif
*/
struct EscSCgrid3 {
    size_t                  n, offs, is, ir, il;
    int                     morder, ma, mb;
    int                     nz, nx, ny, na, nb;
    float                   oz, ox, oy, oa, ob;
    float                   dz, dx, dy, da, db;
    float                   zmin, zmax, xmin, xmax, ymin, ymax, md;
    multi_UBspline_3d_s    *scsplines;
    unsigned char          *mmaped;
    sf_esc_point3           esc_point;
    sf_esc_tracer3          esc_tracer;
    /* Thin plate spline data */
    float                 **L;
    lu_int                 *pvt, ns;
    float                  *bs, *as;
    /* Remote access data for distributed computations */
    bool                    remote;
    int                    *sockets;
    sf_esc_scgrid3_invals  *invals;
    sf_esc_scgrid3_areq    *areqs;
    sf_esc_scgrid3_avals   *avals;
};
/* concrete data type */

/* 8-point stencil, version one:
   b
   ^
   |
   +---+---X---+
   |   |   |   |
   X---X---X---+
   |   |   |   |
   +---X---x---X
   |   |   |   |
   +---X---+---+--->a
*/
/*
#define SCGRID3_TPS_STENCIL1 8
static int sf_scgrid3_tps_ib_stencil1[SCGRID3_TPS_STENCIL1] = 
{ 2, 0, 1, 2, 1, 2, 3, 1 };
static int sf_scgrid3_tps_ia_stencil1[SCGRID3_TPS_STENCIL1] = 
{ 0, 1, 1, 1, 2, 2, 2, 3 };
*/
/* 8-point stencil, version two:
   b
   ^
   |
   X---+---+---X
   |   |   |   |
   +---X---X---+
   |   |   |   |
   +---X---x---+
   |   |   |   |
   X---+---+---X--->a
*/
/*
#define SCGRID3_TPS_STENCIL2 8
static int sf_scgrid3_tps_ib_stencil2[SCGRID3_TPS_STENCIL2] = 
{ 0, 3, 1, 2, 1, 2, 0, 3 };
static int sf_scgrid3_tps_ia_stencil2[SCGRID3_TPS_STENCIL2] = 
{ 0, 0, 1, 1, 2, 2, 3, 3 };
*/
/* 12-point stencil:
   b
   ^
   |
   +---X---X---+
   |   |   |   |
   X---X---X---X
   |   |   |   |
   X---X---x---X
   |   |   |   |
   +---X---X---+--->a
*/

#define SCGRID3_TPS_STENCIL3 12
static int sf_scgrid3_tps_ib_stencil3[SCGRID3_TPS_STENCIL3] = 
{ 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2 };
static int sf_scgrid3_tps_ia_stencil3[SCGRID3_TPS_STENCIL3] = 
{ 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3 };

/* 16-point stencil:
   b
   ^
   |
   X---X---X---X
   |   |   |   |
   X---X---X---X
   |   |   |   |
   X---X---x---X
   |   |   |   |
   X---X---X---X--->a
*/
/*
#define SCGRID3_TPS_STENCIL4 16
static int sf_scgrid3_tps_ib_stencil4[SCGRID3_TPS_STENCIL4] = 
{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
static int sf_scgrid3_tps_ia_stencil4[SCGRID3_TPS_STENCIL4] = 
{ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 };
*/
#define SCGRID3_TPS_MAX_STENCIL 16

/* Initialize the random numbers generator */
static void sf_esc_scgrid3_init_rand (bool devrandom) {
    unsigned int seed;
    struct timeval tv;
    FILE *drand;

    if (false == devrandom || (drand = fopen ("/dev/random","r")) == NULL) {
        gettimeofday (&tv, 0);
        seed = tv.tv_sec + tv.tv_usec;
    } else {
        fread (&seed, sizeof(seed), 1, drand);
        fclose (drand);
    }
    srand (seed);
}

/* Initialize thin-plane spline structures */
static void sf_esc_scgrid3_init_tps (sf_esc_scgrid3 esc_scgrid) {
    const int ns = SCGRID3_TPS_STENCIL3;
    int *bs = sf_scgrid3_tps_ib_stencil3;
    int *as = sf_scgrid3_tps_ia_stencil3;
    lu_int n = ns + 3, info;
    int i;

    esc_scgrid->L = sf_floatalloc2 (ns + 3, ns + 3);
    esc_scgrid->as = sf_floatalloc (ns);
    esc_scgrid->bs = sf_floatalloc (ns);
    esc_scgrid->pvt = (lu_int*)sf_alloc (ns + 3, sizeof(lu_int));

    for (i = 0; i < ns; i++) {
        esc_scgrid->bs[i] = (float)bs[i];
        esc_scgrid->as[i] = (float)as[i];
    }
    for (i = 0; i < (ns + 3); i++) {
        esc_scgrid->pvt[i] = 0;
    }
    sf_tps_build_matrix (esc_scgrid->L, esc_scgrid->bs, esc_scgrid->as, ns, 0.0);
#ifdef MKL_LU
    sgetrf_ (&n, &n, &esc_scgrid->L[0][0], &n, esc_scgrid->pvt, &info);                    
#else
    info = sf_ludlt_decomposition (&esc_scgrid->L[0][0], esc_scgrid->pvt, n);
#endif
    if (info != 0)
        sf_error ("sf_esc_scgrid3_init_tps: TPS matrix inverse failed"); 
    esc_scgrid->ns = ns;
}

sf_esc_scgrid3 sf_esc_scgrid3_init (sf_file scgrid, sf_file scdaemon, sf_esc_tracer3 esc_tracer,
                                    int morder, float frac, bool verb)
/*< Initialize object >*/
{
    size_t nc, nnc = 0;
    int ia, ib, i, j, jj, nab, iab0, iab1, nd, is, on = 1;
    int *nsck, **sck;
    FILE *stream;
    struct sockaddr_in *serv_addr;
    fd_set sset; 
    struct timeval timeout; 
    int valopt, rc, bsiz; 
    socklen_t lon; 

    sf_esc_scgrid3 esc_scgrid = (sf_esc_scgrid3)sf_alloc (1, sizeof (struct EscSCgrid3));

    if (sf_gettype (scgrid) != SF_UCHAR)
        sf_error ("Wrong data type in supercell file");

    if (!sf_histlargeint (scgrid, "n1", (off_t*)&esc_scgrid->n)) sf_error ("No n1= in supercell file");

    if (!sf_histint (scgrid, "Nz", &esc_scgrid->nz)) sf_error ("No Nz= in supercell file");
    if (!sf_histfloat (scgrid, "Dz", &esc_scgrid->dz)) sf_error ("No Dz= in supercell file");
    if (!sf_histfloat (scgrid, "Oz", &esc_scgrid->oz)) sf_error ("No Oz= in supercell file");
    if (!sf_histint (scgrid, "Nx", &esc_scgrid->nx)) sf_error ("No Nx= in supercell file");
    if (!sf_histfloat (scgrid, "Dx", &esc_scgrid->dx)) sf_error ("No Dx= in supercell file");
    if (!sf_histfloat (scgrid, "Ox", &esc_scgrid->ox)) sf_error ("No Ox= in supercell file");
    if (!sf_histint (scgrid, "Ny", &esc_scgrid->ny)) sf_error ("No Ny= in supercell file");
    if (!sf_histfloat (scgrid, "Dy", &esc_scgrid->dy)) sf_error ("No Dy= in supercell file");
    if (!sf_histfloat (scgrid, "Oy", &esc_scgrid->oy)) sf_error ("No Oy= in supercell file");
    if (!sf_histint (scgrid, "Nb", &esc_scgrid->nb)) sf_error ("No Nb= in supercell file");
    if (!sf_histfloat (scgrid, "Db", &esc_scgrid->db)) sf_error ("No Db= in supercell file");
    if (!sf_histfloat (scgrid, "Ob", &esc_scgrid->ob)) sf_error ("No Ob= in supercell file");
    if (!sf_histint (scgrid, "Na", &esc_scgrid->na)) sf_error ("No Na= in supercell file");
    if (!sf_histfloat (scgrid, "Da", &esc_scgrid->da)) sf_error ("No Da= in supercell file");
    if (!sf_histfloat (scgrid, "Oa", &esc_scgrid->oa)) sf_error ("No Oa= in supercell file");

    if (!sf_histfloat (scgrid, "Mdist", &esc_scgrid->md)) esc_scgrid->md = esc_scgrid->dz;

    esc_scgrid->zmin = esc_scgrid->oz;
    esc_scgrid->zmax = esc_scgrid->oz + (esc_scgrid->nz - 1)*esc_scgrid->dz;
    esc_scgrid->xmin = esc_scgrid->ox;
    esc_scgrid->xmax = esc_scgrid->ox + (esc_scgrid->nx - 1)*esc_scgrid->dx;
    esc_scgrid->ymin = esc_scgrid->oy;
    esc_scgrid->ymax = esc_scgrid->oy + (esc_scgrid->ny - 1)*esc_scgrid->dy;

    stream = sf_filestream (scgrid);
    if (stdin == stream)
        sf_error ("Can not mmap supercell file from stdin");

    esc_scgrid->offs = ftello (stream);

    esc_scgrid->scsplines = (multi_UBspline_3d_s*)sf_alloc ((size_t)esc_scgrid->na*(size_t)esc_scgrid->nb,
                                                            sizeof(multi_UBspline_3d_s));
/*
#ifdef DEBUG

    esc_scgrid->mmaped = sf_ucharalloc ((size_t)esc_scgrid->offs +
                                        (size_t)esc_scgrid->n*
                                        (size_t)esc_scgrid->na*
                                        (size_t)esc_scgrid->nb);
    sf_ucharread (esc_scgrid->mmaped, (size_t)esc_scgrid->offs +
                                      (size_t)esc_scgrid->n*
                                      (size_t)esc_scgrid->na*
                                      (size_t)esc_scgrid->nb, scgrid);

#else
*/
    esc_scgrid->mmaped = (unsigned char*)mmap (NULL, (size_t)esc_scgrid->offs +
                                                     (size_t)esc_scgrid->n*
                                                     (size_t)esc_scgrid->na*
                                                     (size_t)esc_scgrid->nb,
                                               PROT_READ, MAP_SHARED,
                                               fileno (stream), 0);
/*
#endif
*/
    nc = esc_scgrid->offs;
    for (ia = 0; ia < esc_scgrid->na; ia++) {
        for (ib = 0; ib < esc_scgrid->nb; ib++) {
            /* Copy spline structure to a separate place to
               make it modifiable */
            memcpy ((void*)&esc_scgrid->scsplines[ia*esc_scgrid->nb + ib],
                    (void*)&esc_scgrid->mmaped[nc], sizeof(multi_UBspline_3d_s));
            nc += sizeof(multi_UBspline_3d_s);
            /* Initialize pointer to coefficients with a correct
               new address in the mmaped area */
            esc_scgrid->scsplines[ia*esc_scgrid->nb + ib].coefs = (float*)&esc_scgrid->mmaped[nc];
            nc += esc_scgrid->scsplines[ia*esc_scgrid->nb + ib].nc;
            nnc += esc_scgrid->scsplines[ia*esc_scgrid->nb + ib].nc;
        }
    }

    if (verb) {
        sf_warning ("Spatial domain dimensions: nz=%d, z=[%g, %g]", esc_scgrid->nz,
                    esc_scgrid->zmin, esc_scgrid->zmax);
        sf_warning ("Spatial domain dimensions: nx=%d, x=[%g, %g]", esc_scgrid->nx,
                    esc_scgrid->xmin, esc_scgrid->xmax);
        sf_warning ("Spatial domain dimensions: ny=%d, y=[%g, %g]", esc_scgrid->ny,
                    esc_scgrid->ymin, esc_scgrid->ymax);
        sf_warning ("Angular domain dimensions: nb=%d, b=[%g, %g]", esc_scgrid->nb,
                    esc_scgrid->ob, esc_scgrid->ob + (esc_scgrid->nb - 1)*esc_scgrid->db);
        sf_warning ("Angular domain dimensions: na=%d, a=[%g, %g]", esc_scgrid->na,
                    esc_scgrid->oa, esc_scgrid->oa + (esc_scgrid->na - 1)*esc_scgrid->da);
        sf_warning ("%lu supercells read, %g Mb of spline coefficients accumulated",
                    (size_t)esc_scgrid->na*(size_t)esc_scgrid->nb, nnc*1e-6);
    }

    esc_scgrid->da *= SF_PI/180.0;
    esc_scgrid->oa *= SF_PI/180.0;
    esc_scgrid->db *= SF_PI/180.0;
    esc_scgrid->ob *= SF_PI/180.0;

    esc_scgrid->is = 0; /* Total number of interpolation steps */
    esc_scgrid->ir = 0; /* Total number of processed points */
    esc_scgrid->il = 0; /* Total number of locally processed points */

    esc_scgrid->morder = morder; /* Interpolation accuracy */

    esc_scgrid->esc_tracer = esc_tracer;
    sf_esc_tracer3_set_mdist (esc_tracer, esc_scgrid->md);
    esc_scgrid->esc_point = sf_esc_point3_init ();

    if (2 != morder) {
        esc_scgrid->L = NULL;
        esc_scgrid->as = NULL;
        esc_scgrid->bs = NULL;
        esc_scgrid->pvt = NULL;
        switch (morder) {
            /* Interpolation stencil length */
            case 3:
                /* Spline interpolation of local 6x6 patch */
                esc_scgrid->ns = 36;
                break;
            case 1:
                /* Bilinear interpolation, 2x2 patch */
                esc_scgrid->ns = 4;
                break;
            case 0:
                /* Nearest neighbor interpolation */
                esc_scgrid->ns = 1;
                break;
        }
    } else {
        /* Thin plate spline; can work with different stencils */
        sf_esc_scgrid3_init_tps (esc_scgrid);
    }

    sf_esc_scgrid3_init_rand (false);

    nc = 0;
    esc_scgrid->sockets = NULL;
    esc_scgrid->ma = 1;
    esc_scgrid->mb = 1;
    if (scdaemon) {
        nab = esc_scgrid->na*esc_scgrid->nb;
        if (!sf_histint (scdaemon, "n2", &nd)) sf_error ("No n2= in supercell daemon file");
        if (!sf_histint (scdaemon, "Ma", &esc_scgrid->ma)) sf_error ("No Ma= in supercell daemon file");
        if (!sf_histint (scdaemon, "Mb", &esc_scgrid->mb)) sf_error ("No Mb= in supercell daemon file");
        serv_addr = sf_alloc (nd, sizeof(struct sockaddr_in));
        sck = (int**)sf_alloc (nab, sizeof(int*)); /* All daemons for an angle */
        nsck = (int*)sf_alloc (nab, sizeof(int)); /* Number of sockets for an angle */
        esc_scgrid->sockets = (int*)sf_alloc (nab, sizeof(int)); /* Connection socket for each angle */
        for (i = 0; i < nab; i++) {
            nsck[i] = 0;
            sck[i] = (int*)sf_alloc (nd, sizeof(int));
        }
        /* Find out how many daemons cover each angle */
        for (j = 0; j < nd; j++) {
            sf_ucharread ((unsigned char*)&serv_addr[j], sizeof (serv_addr[j]), scdaemon);
            sf_ucharread ((unsigned char*)&iab0, sizeof (int), scdaemon);
            sf_ucharread ((unsigned char*)&iab1, sizeof (int), scdaemon);
            if (iab0 > iab1)
                continue;
            /* Store daemon index for this angle patch */
            for (i = iab0; i <= iab1; i++) {
                jj = i;
                if (jj < 0)
                    jj += nab;
                else if (jj >= nab)
                    jj -= nab;
                sck[jj][nsck[jj]] = j;
                nsck[jj]++;
            }
        }
        jj = -1; /* Current daemon index */
        is = -1; /* Current socket */
        for (i = 0; i < nab; i++) {
            j = sck[i][(int)(frac*(float)nsck[i])];
            if (jj == j) { /* Same daemon */
                esc_scgrid->sockets[i] = is;
                continue;
            }
            /* Next daemon, try to establish a connection */
            jj = j;
            /* Create TCP socket */
            if ((is = socket (AF_INET, SOCK_STREAM, 0)) < 0) {
                sf_warning ("socket() failed");
                is = -1;
                continue;
            }
            /* Allow socket descriptor to be reuseable */
            if (setsockopt (is, SOL_SOCKET, SO_REUSEADDR, (char *)&on, sizeof(on)) < 0) {
                sf_warning ("setsockopt() failed");
                close (is);
                is = -1;
                continue;
            }
            on = 1;
            /* Disable TCP buffering of outgoing packets */
#ifdef SOL_TCP
            if (setsockopt (is, SOL_TCP, TCP_NODELAY, (char *)&on, sizeof(on)) < 0) {
#else
            if (setsockopt (is, IPPROTO_TCP, TCP_NODELAY, (char *)&on, sizeof(on)) < 0) {
#endif
                sf_warning ("setsockopt()[TCP_NODELAY] failed");
                close (is);
                is = -1;
                continue;
            }
#ifdef SO_NOSIGPIPE
            on = 1;
            if (setsockopt (new_sd, SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof(on)) < 0 ||) {
                sf_warning ("setsockopt()[SO_NOSIGPIPE] failed");
                close (is);
                continue;
            }
#endif
            /* Set send and receive buffers */
            bsiz = sizeof(sf_esc_scgrid3_areq)*esc_scgrid->ma*esc_scgrid->mb*esc_scgrid->ns;
            if (setsockopt (is, SOL_SOCKET, SO_RCVBUF, &bsiz, sizeof(int)) < 0) {
                sf_warning ("setsockopt()[SO_RCVBUF] failed");
                close (is);
                is = -1;
                continue;
            }
            if (setsockopt (is, SOL_SOCKET, SO_SNDBUF, &bsiz, sizeof(int)) < 0) {
                sf_warning ("setsockopt()[SO_SNDBUF] failed");
                close (is);
                is = -1;
                continue;
            }
            /* Set socket to be non-blocking */
            if (ioctl (is, FIONBIO, (char *)&on) < 0) {
                sf_warning ("ioctl() failed");
                close (is);
                is = -1;
                continue;
            }
            /* Try to establish a connection */
            rc = connect (is, (struct sockaddr *)&serv_addr[j], sizeof(serv_addr[j]));
            if (rc < 0) { 
                if (EINPROGRESS == errno) { 
                    do {
                        /* Connection timeout */
                        timeout.tv_sec = 5*60; 
                        timeout.tv_usec = 0; 
                        FD_ZERO(&sset); 
                        FD_SET(is, &sset); 
                        rc = select (is + 1, NULL, &sset, NULL, &timeout); 
                        if (rc < 0 && errno != EINTR) { 
                            sf_warning ("connect() failed");
                            close (is);
                            is = -1;
                            break; 
                        } else if (rc > 0) { 
                            lon = sizeof(int); 
                            if (getsockopt (is, SOL_SOCKET, SO_ERROR, (void*)(&valopt), &lon) < 0) { 
                                sf_warning ("getsockopt() failed");
                                close (is);
                                is = -1;
                                break;
                            }
                            if (valopt) {
                                sf_warning ("Error in establishing connection for id=%d", j);
                                close (is);
                                is = -1;
                                break;
                            }
                            /* Sucessful connect */
                            break; 
                        } else { 
                            sf_warning ("Connection timeout for id=%d", j);
                            close (is);
                            is = -1;
                            break;
                        } 
                    } while (true);
                } else { 
                    sf_warning ("Connection error for id=%d", j);
                    close (is);
                    is = -1;
                }
            }
            if (is < 0)
                continue;
            esc_scgrid->sockets[i] = is;
            nc++;
        }
        free (serv_addr);
        for (i = 0; i < nab; i++)
            free (sck[i]);
        free (sck);
        free (nsck);
    }
    esc_scgrid->remote = nc != 0;

    esc_scgrid->invals = sf_alloc ((size_t)esc_scgrid->ma*(size_t)esc_scgrid->mb,
                                   sizeof(sf_esc_scgrid3_invals));
    esc_scgrid->areqs = sf_alloc ((size_t)esc_scgrid->ma*(size_t)esc_scgrid->mb*(size_t)esc_scgrid->ns,
                                  sizeof(sf_esc_scgrid3_areq));
    esc_scgrid->avals = sf_alloc ((size_t)esc_scgrid->ma*(size_t)esc_scgrid->mb*(size_t)esc_scgrid->ns,
                                  sizeof(sf_esc_scgrid3_avals));
    if (nc)
        sf_warning ("Using %d remote daemons for accessing escape solutions", nc);
    else
        sf_warning ("Using local data for accessing escape solutions");

    return esc_scgrid;
}

/* Disconnect from a socket */
static void sf_cram_scgrid3_disconnect (sf_esc_scgrid3 esc_scgrid, int is) {
    int i;

    close (is);

    /* Remove socket id from all angle references */
    for (i = 0; i < esc_scgrid->na*esc_scgrid->nb; i++) {
        if (esc_scgrid->sockets[i] == is)
            esc_scgrid->sockets[i] = -1;
    }
}

void sf_esc_scgrid3_close (sf_esc_scgrid3 esc_scgrid, bool verb)
/*< Destroy object >*/
{
    int iab, nab = esc_scgrid->na*esc_scgrid->nb;
    if (verb)
        sf_warning ("%lu points processed, %g interpolation steps per point performed",
                    esc_scgrid->ir, (float)esc_scgrid->is/(float)esc_scgrid->ir);
        sf_warning ("%lu points processed locally (%g%%)",
                    esc_scgrid->il, 100.0*(float)esc_scgrid->il/(float)esc_scgrid->ir);
    /* Close all existing connections */
    if (esc_scgrid->sockets) {
        for (iab = 0; iab < nab; iab++) {
            if (esc_scgrid->sockets[iab] != -1)
                sf_cram_scgrid3_disconnect (esc_scgrid, esc_scgrid->sockets[iab]);
        }
        free (esc_scgrid->sockets);
    }
    if (2 == esc_scgrid->morder) {
        /* Release special structures used for thin plate spline interpolation */
        free (esc_scgrid->L[0]);
        free (esc_scgrid->L);
        free (esc_scgrid->as);
        free (esc_scgrid->bs);
        free (esc_scgrid->pvt);
    }
    free (esc_scgrid->invals);
    free (esc_scgrid->areqs);
    free (esc_scgrid->avals);
/*
#ifdef DEBUG
    free (esc_scgrid->mmaped);
#else
*/
    munmap (esc_scgrid->mmaped, (size_t)esc_scgrid->offs +
                                (size_t)esc_scgrid->n*
                                (size_t)esc_scgrid->na*
                                (size_t)esc_scgrid->nb);
/*
#endif
*/
    free (esc_scgrid->scsplines);
    sf_esc_point3_close (esc_scgrid->esc_point);
    free (esc_scgrid);
}

/* Compute one value locally */
static void sf_cram_scgrid3_get_lvalue (sf_esc_scgrid3 esc_scgrid, sf_esc_scgrid3_areq *areq,
                                        sf_esc_scgrid3_avals *aval) {
    eval_multi_UBspline_3d_s (&esc_scgrid->scsplines[areq->iab],
                              areq->y, areq->x, areq->z, aval->vals);
    aval->ud1 = areq->ud1;
    aval->ud2 = areq->ud2;
}

/* Compute multiple values either locally or remotely */
static void sf_cram_scgrid3_get_values (sf_esc_scgrid3 esc_scgrid, sf_esc_scgrid3_areq *areqs,
                                        sf_esc_scgrid3_avals *avals, int n) {
    int i, ii, ie, is, mis = -1, iab;
    int len = 0, rc, ns = 0, desc_ready;
    size_t id;
    fd_set sset, wset;
    struct timeval timeout;

    if (0 == n)
        return;

    /* Extract values locally */
    if (false == esc_scgrid->remote) {
        for (i = 0; i < n; i++) {
            sf_cram_scgrid3_get_lvalue (esc_scgrid, &areqs[i], &avals[i]);
            esc_scgrid->il++;
        }
        return;
    }

    /* Message ID for this batch of points */
    id = 0;
    while (0 == id && id >= (size_t)(-n))
        id = rand ()*rand ();
    FD_ZERO(&sset);

    /* Skip patches without receivers */
    ii = 0;
    while (ii < n && -1 == esc_scgrid->sockets[areqs[ii].iab]) {
        areqs[ii].id = id;
        ii++;
    }
    if (ii == n)
        sf_error ("Lost all connections");

    /* Socket */
    is = esc_scgrid->sockets[areqs[ii].iab];
    if (is > mis)
        mis = is;
    ie = ii; /* Last and first requests in the patch */
    iab = areqs[ie].iab;
    /* Find maximum consecutive patch of angles to
       be sent to a single socket */
    do {
        while (ie < n && (iab == areqs[ie].iab ||
                          is == esc_scgrid->sockets[areqs[ie].iab])) {
            iab = areqs[ie].iab;
            areqs[ie].id = id + ie;
            ie++;
        }
        /* Send loop for this patch */
        len = 0;
        while (len < sizeof(sf_esc_scgrid3_areq)*(ie - ii)) {
            rc = send (is, (const void*)(((unsigned char*)&areqs[ii]) + len),
                       sizeof(sf_esc_scgrid3_areq)*(ie - ii) - len, MSG_NOSIGNAL);
            if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                sf_warning ("Can not send data for iab=[%d - %d], disconnecting",
                            areqs[ii].iab, areqs[ie].iab - 1);
                sf_cram_scgrid3_disconnect (esc_scgrid, is);
                len = 0;
                break;
            }
            if (rc > 0)
                len += rc;
        }
        if (len && !FD_ISSET (is, &sset))
            FD_SET(is, &sset);
        ns += len/sizeof(sf_esc_scgrid3_areq);
        if (ie < n) { /* Next patch */
            ii = ie;
            /* Find next socket */
            while (ii < n && -1 == esc_scgrid->sockets[areqs[ii].iab]) {
                areqs[ii].id = id; /* Skip requests without receivers */
                ii++;
            }
            ie = ii;
            if (ie < n) {
                is = esc_scgrid->sockets[areqs[ii].iab];
                if (is > mis)
                    is = mis;
                iab = areqs[ie].iab;
            }
        }
    } while (ie < n); /* Done sending */
    /* Now, try to get results back */
    ii = 0;
    while (ns) { /* Poll sockets for incoming data */
        /* Wait for 60 secs max */
        timeout.tv_sec  = 60;
        timeout.tv_usec = 0;
        memcpy (&wset, &sset, sizeof(sset));
        rc = select (mis + 1, &wset, NULL, NULL, &timeout);
        if (0 == rc)
            break;
        if (rc < 0)
            sf_error ("select() failed");
        desc_ready = rc;
        for (is = 0; is <= mis && desc_ready > 0; is++) {
            /* Check to see if this descriptor is ready */
            if (!FD_ISSET (is, &wset))
                continue;
            desc_ready--;
            len = 0;
            /* Receive loop */
            while (len < sizeof(sf_esc_scgrid3_avals)*ns) {
                /* Receive job result from the client */
                rc = recv (is, (void*)(((unsigned char*)&avals[ii]) + len),
                           sizeof(sf_esc_scgrid3_avals)*ns - len, 0);
                if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                    /* Connection terminated */
                    if (0 == rc)
                        sf_warning ("The server has closed connection for socket %d", is);
                    else
                        sf_warning ("Can not receive data for socket %d, disconnecting", is);
                    sf_cram_scgrid3_disconnect (esc_scgrid, is);
                    FD_CLR(is, &sset);
                    len = 0;
                    break;
                }
                if (rc > 0)
                    len += rc;
                if (0 == len % sizeof(sf_esc_scgrid3_avals)) /* A few complete patches have been received */
                    break;
            }
            if (len > 0) {
                if (len % sizeof(sf_esc_scgrid3_avals)) {
                    sf_warning ("Partial receive from socket %d\n", i);
                    len = 0;
                } else {
                    len /= sizeof(sf_esc_scgrid3_avals);
                    for (i = 0; i < len; i++) {
                        if (avals[ii].id >= id && avals[ii].id < (id + (size_t)n)) {
                            if (avals[ii].vals[0] != SF_HUGE) {
                                areqs[avals[ii].id - id].id = 0; /* Mark as received */
                                ns--;
                                ii++;
                            } else {
                                sf_warning ("Server replied that the angle is out of bounds");
                                break;
                            }
                        } else {
                            sf_warning ("Received garbage from socket %d", is);
                            break;
                        }
                    }
                }
            }
        } /* Loop over ready sockets */
    } /* End of polling */
    if (ns)
        sf_warning ("Timeout, %d angles are left", ns);

    /* Extract values locally */
    for (i = 0; i < n; i++) {
        if (areqs[i].id != 0) {
            sf_cram_scgrid3_get_lvalue (esc_scgrid, &areqs[i], &avals[ii]);
            esc_scgrid->il++;
            ii++;
        }
    }
}

/* Compute interger index for a float value v according to sampling df
   and zero at f0, store fractional part in f */
static int sf_cram_scgrid3_ftoi (float v, float f0, float df, float *f) {
    int i;
    float ff;

    ff = (v - f0)/df;
    i = floorf (ff);
    if (f)
        *f = ff - (float)i;
    return i;
}

/* Convert normalized phase vector components into polar angles */
static void sf_esc_scgrid3_p_to_ab (float pz, float px, float py,
                                    float *a, float *b) {
    float l;

    l = sqrt (px*px + py*py);
    if (l > 1e-6) {
        if (py >= 0.0)
            *a = acosf (px/l);
        else
            *a = 2.0*SF_PI - acosf (px/l);
    }
    l = sqrt (pz*pz + px*px + py*py);
    *b = acosf (pz/l);
}

/* Convert quaternion of rotation into a new phase direction */
static void sf_esc_scgrid3_q_to_ab (float *q, float *a, float *b) {
    int i, j;
    float vf[3], vt[3];
    float M[3][3];

    vf[0] = cosf (*b);
    vf[1] = sinf (*b)*cosf (*a);
    vf[2] = sinf (*b)*sinf (*a);

    /* Obtain rotation to the new phase position */
    sf_quat_norm (q);
    sf_quat_rotmat (q, &M[0][0]);

    /* Rotate the old phase vector */
    for (i = 0; i < 3; i++) {
        vt[i] = 0.0;
        for (j = 0; j < 3; j++) {
            vt[i] += M[i][j]*vf[j];
        }
    }
    /* Convert to azimuth and inclination */
    sf_esc_scgrid3_p_to_ab (vt[0], vt[1], vt[2], a, b);
}

/* Flip azimuth by +/-PI (this happens when inclination is moved back to
   its interval of [0;PI) */
static int sf_esc_scgrid3_flip_ia (sf_esc_scgrid3 esc_scgrid, int ia) {
    if (ia < esc_scgrid->na/2)
        return ia + esc_scgrid->na/2;
    else
        return ia - esc_scgrid->na/2;
}

/* Adjust azimuth according to its periodicity: [0;2PI) for the azimuth */
static int sf_esc_scgrid3_bound_ia (sf_esc_scgrid3 esc_scgrid, int ia) {
    if (ia < 0) {
        ia += esc_scgrid->na;
    } else if (ia >= esc_scgrid->na) {
        ia -= esc_scgrid->na;
    }
    return ia;
}

/* Adjust azimuth and inclination angles according to their periodicity:
   [0;PI) for the inclination, [0;2PI) for the azimuth */
static void sf_esc_scgrid3_bound_iaib (sf_esc_scgrid3 esc_scgrid, int *ia, int *ib) {
    /* Azimuth */
    *ia = sf_esc_scgrid3_bound_ia (esc_scgrid, *ia);
    /* Inclination */
    if (*ib < 0) {
        *ib = -(*ib + 1);
        *ia = sf_esc_scgrid3_flip_ia (esc_scgrid, *ia);
    }
    if (*ib >= esc_scgrid->nb) {
        *ib = 2*esc_scgrid->nb - (*ib) - 1;
        *ia = sf_esc_scgrid3_flip_ia (esc_scgrid, *ia);
    }
}

/* Perform nearest neighbor interpolation of the local escape solution
   from stencil values (avals), update the previous iteration (inval) */
static void sf_esc_scgrid3_nearest_interp (sf_esc_scgrid3 esc_scgrid,
                                           sf_esc_scgrid3_invals *inval,
                                           sf_esc_scgrid3_avals *avals) {
    int i;

    for (i = 0; i < ESC3_NUM; i++)
        inval->vals[i] += avals->vals[i];
    sf_esc_scgrid3_q_to_ab (&avals->vals[ESC3_NUM], &inval->a, &inval->b);
}

/* Perform bilinear interpolation of the local escape solution
   from stencil values (avals), update the previous iteration (inval) */
static void sf_esc_scgrid3_bilinear_interp (sf_esc_scgrid3 esc_scgrid,
                                          sf_esc_scgrid3_invals *inval,
                                          sf_esc_scgrid3_avals *avals) {
    int i;
    float v[ESC3_NUM + 3];

    /* Bilinear interpolation */
    for (i = 0; i < (ESC3_NUM + 4); i++) {
        v[i] = avals[0].vals[i]*(1.0 - inval->fb)*(1.0 - inval->fa) +
               avals[1].vals[i]*inval->fb*(1.0 - inval->fa) +
               avals[2].vals[i]*(1.0 - inval->fb)*inval->fa +
               avals[3].vals[i]*inval->fb*inval->fa;
    }
    /* Update previous values */
    for (i = 0; i < ESC3_NUM; i++)
        inval->vals[i] += v[i];
    sf_esc_scgrid3_q_to_ab (&v[ESC3_NUM], &inval->a, &inval->b);
}

/* Perform spline interpolation of the local escape solution
   from stencil values (avals), update the previous iteration (inval) */
static void sf_esc_scgrid3_spline_interp (sf_esc_scgrid3 esc_scgrid,
                                          sf_esc_scgrid3_invals *inval,
                                          sf_esc_scgrid3_avals *avals) {
    int i, k;
    float f[ESC3_NUM + 4][36], v[ESC3_NUM + 3];
    Ugrid z_grid, x_grid;
    BCtype_s zBC, xBC;
    multi_UBspline_2d_s *escspline = NULL;

    /* Put escape functions into separate arrays */ 
    for (k = 0; k < (ESC3_NUM + 4); k++) {
        for (i = 0; i < 36; i++) {
            f[k][i] = avals[i].vals[k];
        }
    }
    /* Spline interpolation */
    z_grid.start = 0; z_grid.end = 5; z_grid.num = 6;
    x_grid.start = 0; x_grid.end = 5; x_grid.num = 6;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    escspline = create_multi_UBspline_2d_s (x_grid, z_grid, xBC, zBC, ESC3_NUM + 4);
    for (k = 0; k < (ESC3_NUM + 4); k++) {
        set_multi_UBspline_2d_s (escspline, k, &f[k][0]);
    }
    eval_multi_UBspline_2d_s (escspline, 2.0 + inval->fa, 2.0 + inval->fb, v);
    /* Update previous values */
    for (i = 0; i < ESC3_NUM; i++)
        inval->vals[i] += v[i];
    sf_esc_scgrid3_q_to_ab (&v[ESC3_NUM], &inval->a, &inval->b);
    destroy_Bspline (escspline);
}

/* Perform thin plate spline interpolation of the local escape solution
   from stencil values (avals), update the previous iteration (inval) */
static void sf_esc_scgrid3_tps_interp (sf_esc_scgrid3 esc_scgrid,
                                       sf_esc_scgrid3_invals *inval,
                                       sf_esc_scgrid3_avals *avals) {
    int i, j;
    float f[SCGRID3_TPS_MAX_STENCIL + 3], w[SCGRID3_TPS_MAX_STENCIL + 3],
          v[ESC3_NUM + 4];
    bool overwrt = false;
    lu_int n = esc_scgrid->ns + 3, info;
#ifdef MKL_LU
    lu_int nrhs = 1, ldb = n
#endif

    /* Loop over escape functions */
    for (i = 0; i < (ESC3_NUM + 4); i++) {
        /* Put one escape function type into a contiguous array */ 
        for (j = 0; j < esc_scgrid->ns; j++) {
            f[j] = avals[j].vals[i];
        }
        for (j = esc_scgrid->ns; j < (esc_scgrid->ns + 3); j++) {
            f[j] = 0.0;
        }
        /* Compute thin-plate spline coefficients */
#ifdef MKL_LU
        overwrt = true;
        sgetrs_ ("N", &n, &nrhs, &esc_scgrid->L[0][0], &n, esc_scgrid->pvt, f, &ldb, &info);
#else
        info = sf_ludtl_solve (&esc_scgrid->L[0][0], f, esc_scgrid->pvt, w, n);
#endif
        if (info != 0)
            sf_error ("sf_esc_scgrid3_tps_interp: linear system solution error");
        /* Do interpolation */
        v[i] = sf_tps_compute_point (overwrt ? f : w, esc_scgrid->bs, esc_scgrid->as,
                                     esc_scgrid->ns, 1.0 + inval->fb, 1.0 + inval->fa);
        /* Update previous values */
        if (i < ESC3_NUM)
            inval->vals[i] += v[i];
    }
    sf_esc_scgrid3_q_to_ab (&v[ESC3_NUM], &inval->a, &inval->b);
}

/* Perform interpolation of the local escape solution from stencil values (avals)
   according to the current interpolation order, update the previous iteration (inval) */
static void sf_esc_scgrid3_interp_point (sf_esc_scgrid3 esc_scgrid,
                                         sf_esc_scgrid3_invals *inval,
                                         sf_esc_scgrid3_avals *avals) {
    if (3 == esc_scgrid->morder)
        sf_esc_scgrid3_spline_interp (esc_scgrid, inval, avals);
    else if (2 == esc_scgrid->morder)
        sf_esc_scgrid3_tps_interp (esc_scgrid, inval, avals);
    else if (0 == esc_scgrid->morder)
        sf_esc_scgrid3_nearest_interp (esc_scgrid, inval, avals);
    else
        sf_esc_scgrid3_bilinear_interp (esc_scgrid, inval, avals);
}

/* Prepare phase space values request (areq) accoring to the point (aval)
   for the nearest interpolation */
static void sf_esc_scgrid3_nearest_prepare (sf_esc_scgrid3 esc_scgrid, int iia, int iib,
                                            sf_esc_scgrid3_invals *aval,
                                            sf_esc_scgrid3_areq *areq) {
    int ia, ib;

    ia = iia;
    ib = iib;
    sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
    areq->iab = ia*esc_scgrid->nb + ib;
}

/* Prepare phase space values request (areq) accoring to the point (aval)
   for the spline interpolation stencil */
static void sf_esc_scgrid3_bilinear_prepare (sf_esc_scgrid3 esc_scgrid, int iia, int iib,
                                             sf_esc_scgrid3_invals *aval,
                                             sf_esc_scgrid3_areq *areq) {
    int i, j, ia, ib;

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            ia = iia + i;
            ib = iib + j;
            sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
            areq[i].iab = ia*esc_scgrid->nb + ib;
        }
    }
}

/* Prepare phase space values request (areq) accoring to the point (aval)
   for the spline interpolation stencil */
static void sf_esc_scgrid3_spline_prepare (sf_esc_scgrid3 esc_scgrid, int iia, int iib,
                                           sf_esc_scgrid3_invals *aval,
                                           sf_esc_scgrid3_areq *areq) {
    int i, j, ia, ib;

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            ia = iia + i - 2;
            ib = iib + j - 2;
            sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
            areq[i].iab = ia*esc_scgrid->nb + ib;
        }
    }
}

/* Prepare phase space values request (areq) accoring to the point (aval)
   for the thin plate spline interpolation stencil */
static void sf_esc_scgrid3_tps_prepare (sf_esc_scgrid3 esc_scgrid, int iia, int iib,
                                        sf_esc_scgrid3_invals *aval,
                                        sf_esc_scgrid3_areq *areq) {
    int i, ia, ib;

    for (i = 0; i < esc_scgrid->ns; i++) {
        ib = esc_scgrid->bs[i] - 1 + iib;
        ia = esc_scgrid->as[i] - 1 + iia;
        sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
        areq[i].iab = ia*esc_scgrid->nb + ib;
    }
}


/* Prepare phase space values request (areq) accoring to the point (aval)
   and the current interpolation order */
static void sf_esc_scgrid3_prepare_request (sf_esc_scgrid3 esc_scgrid, int in,
                                            sf_esc_scgrid3_invals *aval,
                                            sf_esc_scgrid3_areq *areq) {
    int i, iib, iia;
    /* ib, ia */
    iib = sf_cram_scgrid3_ftoi (aval->b, esc_scgrid->ob, esc_scgrid->db, &aval->fb);
    iia = sf_cram_scgrid3_ftoi (aval->a, esc_scgrid->oa, esc_scgrid->da, &aval->fa);

    if (3 == esc_scgrid->morder)
        sf_esc_scgrid3_spline_prepare (esc_scgrid, iia, iib, aval, areq);
    else if (2 == esc_scgrid->morder)
        sf_esc_scgrid3_tps_prepare (esc_scgrid, iia, iib, aval, areq);
    else if (0 == esc_scgrid->morder)
        sf_esc_scgrid3_nearest_prepare (esc_scgrid, iia, iib, aval, areq);
    else
        sf_esc_scgrid3_bilinear_prepare (esc_scgrid, iia, iib, aval, areq);

    for (i = 0; i < esc_scgrid->ns; i++) {
        areq[i].z = aval->vals[ESC3_Z];
        areq[i].x = aval->vals[ESC3_X];
        areq[i].y = aval->vals[ESC3_Y];
        areq[i].ud1 = in; /* Index in the outer input array */
        areq[i].ud2 = i; /* Stencil index */
    }

}

#define SPEPS 1e-3

/* Check if the point (avals) is inside the phase space supercell grid,
   if not, then project it to the global boundary */
static bool sf_esc_scgrid3_is_inside (sf_esc_scgrid3 esc_scgrid,
                                      float *aval, float *a, float *b, bool project)
{
    if ((aval[ESC3_Z] <= esc_scgrid->zmin + SPEPS*esc_scgrid->dz) ||
        (aval[ESC3_Z] >= esc_scgrid->zmax - SPEPS*esc_scgrid->dz) ||
        (aval[ESC3_X] <= esc_scgrid->xmin + SPEPS*esc_scgrid->dx) ||
        (aval[ESC3_X] >= esc_scgrid->xmax - SPEPS*esc_scgrid->dx) ||
        (aval[ESC3_Y] <= esc_scgrid->ymin + SPEPS*esc_scgrid->dy) ||
        (aval[ESC3_Y] >= esc_scgrid->ymax - SPEPS*esc_scgrid->dy)) {
        /* Do ray tracing, if the point is outside of the supergrid bounds */
        if (project && sf_esc_tracer3_inside (esc_scgrid->esc_tracer, &aval[ESC3_Z],
                                              &aval[ESC3_X], &aval[ESC3_Y], true)) {
            sf_esc_tracer3_compute (esc_scgrid->esc_tracer, aval[ESC3_Z], aval[ESC3_X],
                                    aval[ESC3_Y], *b, *a, aval[ESC3_T],
#ifdef ESC_EQ_WITH_L
                                    aval[ESC3_L],
#else
                                    0.0,
#endif
                                    esc_scgrid->esc_point, a, b);
            aval[ESC3_Z] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_Z);
            aval[ESC3_X] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_X);
            aval[ESC3_Y] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_Y);
            aval[ESC3_T] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_T);
#ifdef ESC_EQ_WITH_L
            aval[ESC3_L] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_L);
#endif
        }
        return false;
    }
    return true;
}

/* Comparison routine for qsort() by angle indices */
static int sf_esc_scgrid3_areq_sort (const void *v1, const void *v2) {
    sf_esc_scgrid3_areq *areq1 = (sf_esc_scgrid3_areq*)v1;
    sf_esc_scgrid3_areq *areq2 = (sf_esc_scgrid3_areq*)v2;

    if (areq1->iab < areq2->iab)
        return -1;
    else if (areq1->iab > areq2->iab)
        return 1;
    else
        return 0;
}

/* Comparison routine for qsort() of avals by point and stencil indices */
static int sf_esc_scgrid3_avals_sort (const void *v1, const void *v2) {
    sf_esc_scgrid3_avals *aval1 = (sf_esc_scgrid3_avals*)v1;
    sf_esc_scgrid3_avals *aval2 = (sf_esc_scgrid3_avals*)v2;

    /* Point index */
    if (aval1->ud1 < aval2->ud1)
        return -1;
    else if (aval1->ud1 > aval2->ud1)
        return 1;
    else {
        /* Stencil index */
        if (aval1->ud2 < aval2->ud2)
            return -1;
        else if (aval1->ud2 > aval2->ud2)
            return 1;
        else       
            return 0;
    }
}


void sf_esc_scgrid3_compute (sf_esc_scgrid3 esc_scgrid, float z, float x, float y,
                             float a0, float da, float b0, float db, int na, int nb,
                             float *avals)
/*< Compute escape values for a point with subsurface coordinates (z, x, y, b, a)
    by stitching local escape solutions in supercells of a phase-space grid >*/
{
    int i, ii, ie, io, in, ia, ib, iap, ibp;
    int nap = na/esc_scgrid->ma, nbp = nb/esc_scgrid->mb;
    sf_esc_scgrid3_invals *input = esc_scgrid->invals;
    sf_esc_scgrid3_areq *areqs = esc_scgrid->areqs;
    sf_esc_scgrid3_avals *output = esc_scgrid->avals;

    if (na % esc_scgrid->ma)
        sf_warning ("na should be divisible by ma");
    if (nb % esc_scgrid->mb)
        sf_warning ("nb should be divisible by mb");

    /* Check if can process this point */
    avals[ESC3_Z] = z;
    avals[ESC3_X] = x;
    avals[ESC3_Y] = y;
    if (false == sf_esc_scgrid3_is_inside (esc_scgrid, avals, NULL, NULL, false)) {
        /* Just trace it to the boundary then */
        for (ia = 0; ia < na; ia++) {
            for (ib = 0; ib < nb; ib++) {
                sf_esc_tracer3_compute (esc_scgrid->esc_tracer, z, x, y,
                                        b0 + ib*db, a0 + ia*da, 0.0, 0.0,
                                        esc_scgrid->esc_point, NULL, NULL);
                i = (ia*nb + ib)*ESC3_NUM;
                avals[i + ESC3_Z] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_Z);
                avals[i + ESC3_X] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_X);
                avals[i + ESC3_Y] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_Y);
                avals[i + ESC3_T] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_T);
#ifdef ESC_EQ_WITH_L
                avals[i + ESC3_L] = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_L);
#endif                
            }
        }
        return;
    }

    for (iap = 0; iap < nap; iap++) { /* Loop over patches */
        for (ibp = 0; ibp < nbp; ibp++) { /* Loop over patches */
            /* Initialize starting positions in the phase space */
            for (ia = 0; ia < esc_scgrid->ma; ia++) {
                for (ib = 0; ib < esc_scgrid->mb; ib++) {
                    input[ia*esc_scgrid->mb + ib].vals[ESC3_Z] = z;
                    input[ia*esc_scgrid->mb + ib].vals[ESC3_X] = x;
                    input[ia*esc_scgrid->mb + ib].vals[ESC3_Y] = y;
                    input[ia*esc_scgrid->mb + ib].vals[ESC3_T] = 0.0;
#ifdef ESC_EQ_WITH_L
                    input[ia*esc_scgrid->mb + ib].vals[ESC3_L] = 0.0;
#endif
                    input[ia*esc_scgrid->mb + ib].b = b0 + (ibp*esc_scgrid->mb + ib)*db;
                    input[ia*esc_scgrid->mb + ib].a = a0 + (iap*esc_scgrid->ma + ia)*da;
                    /* Output index */
                    input[ia*esc_scgrid->mb + ib].iab = (iap*esc_scgrid->ma + ia)*nb +
                                                        (ibp*esc_scgrid->mb + ib);
                }
            }
            /* Number of elements in the input array */
            in = esc_scgrid->ma*esc_scgrid->mb;
            /* First and last indices in the input array */
            ii = 0;
            ie = esc_scgrid->ma*esc_scgrid->mb - 1;

            do {
                io = 0;
                for (i = ii; i <= ie; i++) {
                    if (input[i].iab >= 0) {
                        sf_esc_scgrid3_prepare_request (esc_scgrid, i, &input[i],
                                                        &areqs[io*esc_scgrid->ns]);
                        io++;
                    }
                }
                /* Sort the outgoing requests by angle index */
                qsort (areqs, in*esc_scgrid->ns, sizeof(sf_esc_scgrid3_areq),
                       sf_esc_scgrid3_areq_sort);

                /* Extract values from across the phase space; remotely or locally */
                sf_cram_scgrid3_get_values (esc_scgrid, areqs, output, in*esc_scgrid->ns);

                /* Sort the incoming angle values by point and stencil index */
                qsort (output, in*esc_scgrid->ns, sizeof(sf_esc_scgrid3_avals),
                       sf_esc_scgrid3_avals_sort);

                /* Find completed points */
                io = in;
                for (i = 0; i < io; i++) {
                    sf_esc_scgrid3_interp_point (esc_scgrid, &input[output[i*esc_scgrid->ns].ud1],
                                                 &output[i*esc_scgrid->ns]);
                    if (false == sf_esc_scgrid3_is_inside (esc_scgrid, input[output[i*esc_scgrid->ns].ud1].vals,
                                                           &input[output[i*esc_scgrid->ns].ud1].a,
                                                           &input[output[i*esc_scgrid->ns].ud1].b, true)) {
                        /* Move completed points to the output array */
                        memcpy (&avals[input[output[i*esc_scgrid->ns].ud1].iab*ESC3_NUM],
                                input[output[i*esc_scgrid->ns].ud1].vals, sizeof(float)*ESC3_NUM);
                        in--;
                        input[output[i*esc_scgrid->ns].ud1].iab = -1;
                        esc_scgrid->ir++;
                    } else
                        esc_scgrid->is++;
                }
                /* Find new first and last indices in the input array */
                while (input[ii].iab == -1 && ii < ie)
                    ii++;
                while (input[ie].iab == -1 && ie > ii)
                    ie--;
            } while (in);
        }
    }        
}

