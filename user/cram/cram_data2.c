/* Common-shot data for angle migration (both 2-D and 3-D) */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/socket.h>
#ifdef sun
#include <sys/filio.h>
#else
#include <sys/ioctl.h>
#endif
#include <sys/select.h>
#include <unistd.h>
#include <netinet/in.h>
#include <netinet/tcp.h>

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

#ifndef SOL_TCP
#define SOL_TCP IPPROTO_TCP
#endif

#include <rsf.h>

#ifndef _cram_data2_h

typedef struct CRAMData2 *sf_cram_data2;
/* abstract data type */
/*^*/

typedef struct {
    size_t id;    /* Message ID */
    size_t i, n;  /* Index of the first trace, number of traces */
} sf_cram_data_trreq;
/* Structure for requesting several traces */
/*^*/

typedef struct {
    size_t id;        /* Message ID */
    size_t n;         /* Number of received traces */
    float *samples;   /* Trace samples */
} sf_cram_data_trvals;
/* Structure for getting requested traces back */
/*^*/

#define CRAM_TRBUF 10
/*^*/

#endif

struct CRAMData2 {
    int                  nt;
    size_t               i, ntr, n, offs, nb;
    float                t0, dt, trd;
    bool                 kmah, filter, erefl;
    float               *data;
    char                *mmaped;
    FILE                *stream;
    sf_cram_data_trvals *trvals;
    size_t               nc;
    int                 *sockets;
};
/* concrete data type */

sf_cram_data2 sf_cram_data2_init (sf_file data, sf_file ddaemon)
/*< Initialize object >*/
{
    size_t n, j, i0, i1;
    int is, on = 1, valopt, rc, bsiz, nd;
    struct sockaddr_in *serv_addr;
    fd_set sset; 
    struct timeval timeout; 
    socklen_t lon; 

    sf_cram_data2 cram_data = (sf_cram_data2)sf_alloc (1, sizeof (struct CRAMData2));

    /* Filter data according to slope */
    if (!sf_histbool (data, "filter", &cram_data->filter)) sf_error ("No filter= in data");
    /* Use KMAH phase shifts */
    if (!sf_histbool (data, "KMAH", &cram_data->kmah)) sf_error ("No KMAH= in data");
    /* Use exploding reflector assumption */
    if (!sf_histbool (data, "ExplRefl", &cram_data->erefl)) cram_data->erefl = false;

    if (!sf_histint (data, "n1", &cram_data->nt)) sf_error ("No n1= in data");
    if (!sf_histfloat (data, "d1", &cram_data->dt)) sf_error ("No d1= in data");
    if (!sf_histfloat (data, "o1", &cram_data->t0)) sf_error ("No o1= in data");

    cram_data->n = sf_leftsize (data, 1);

    /* Total size of the file memory mapping */
    n = (size_t)cram_data->nt*(size_t)cram_data->n*(size_t)sizeof(float);

    if (cram_data->kmah)
        sf_warning ("Using data prepared for KMAH shifts");
    if (cram_data->filter)
        sf_warning ("Using data prepared for anti-aliasing filter");
    if (cram_data->erefl)
        sf_warning ("Assuming data modeled by exploding reflector");
    sf_warning ("Total data size: %g Mb", 1e-6*(float)n);

    if (cram_data->kmah)
        cram_data->n /= (size_t)2;

    cram_data->stream = sf_filestream (data);
    cram_data->offs = ftello (cram_data->stream);

    cram_data->nc = 0;
    cram_data->nb = 0;
    cram_data->sockets = NULL;
    cram_data->trvals = NULL;

    if (ddaemon) { /* Try to establish a remote access mode */
        if (!sf_histint (ddaemon, "n2", &nd)) sf_error ("No n2= in data daemon file");
        if (!sf_histfloat (ddaemon, "trd", &cram_data->trd)) sf_error ("No trd= in data daemon file");
        serv_addr = sf_alloc (nd, sizeof(struct sockaddr_in));
        cram_data->sockets = (int*)sf_alloc (nd, sizeof(int));
        for (j = 0; j < nd; j++) {
            sf_ucharread ((unsigned char*)&serv_addr[j], sizeof(serv_addr[j]), ddaemon);
            sf_ucharread ((unsigned char*)&i0, sizeof(size_t), ddaemon);
            sf_ucharread ((unsigned char*)&i1, sizeof(size_t), ddaemon);
            if (i0 > i1)
                continue;
            /* Create TCP socket */
            if ((is = socket (AF_INET, SOCK_STREAM, 0)) < 0) {
                sf_warning ("socket() failed");
                continue;
            }
            /* Allow socket descriptor to be reuseable */
            if (setsockopt (is, SOL_SOCKET, SO_REUSEADDR, (char *)&on, sizeof(on)) < 0) {
                sf_warning ("setsockopt() failed");
                close (is);
                continue;
            }
            on = 1;
            /* Disable TCP buffering of outgoing packets */
#ifndef TCP_CORK
            if (setsockopt (is, SOL_TCP, TCP_NODELAY, (char *)&on, sizeof(on)) < 0) {
                sf_warning ("setsockopt()[TCP_NODELAY] failed");
                close (is);
                continue;
            }
#endif
#ifdef SO_NOSIGPIPE
            on = 1;
            if (setsockopt (is, SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof(on)) < 0) {
                sf_warning ("setsockopt()[SO_NOSIGPIPE] failed");
                close (is);
                continue;
            }
#endif
            /* Set receive buffer size */
            bsiz = cram_data->kmah ? sizeof(sf_cram_data_trvals) + CRAM_TRBUF*cram_data->nt*sizeof(float)*2
                                   : sizeof(sf_cram_data_trvals) + CRAM_TRBUF*cram_data->nt*sizeof(float);
            if (setsockopt (is, SOL_SOCKET, SO_RCVBUF, &bsiz, sizeof(int)) < 0) {
                sf_warning ("setsockopt()[SO_RCVBUF] failed");
                close (is);
                continue;
            }
            /* Set socket to be non-blocking */
            if (ioctl (is, FIONBIO, (char *)&on) < 0) {
                sf_warning ("ioctl() failed");
                close (is);
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
                            /* Successful connect */
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
            cram_data->sockets[cram_data->nc] = is;
            cram_data->nc++;
        }
        free (serv_addr);
    }

    if (cram_data->nc) {
        if (cram_data->kmah)
            cram_data->trvals = sf_alloc (1, sizeof(sf_cram_data_trvals) +
                                             CRAM_TRBUF*cram_data->nt*sizeof(float)*2);
        else
            cram_data->trvals = sf_alloc (1, sizeof(sf_cram_data_trvals) +
                                             CRAM_TRBUF*cram_data->nt*sizeof(float));
        cram_data->data = &cram_data->trvals->samples[0];
        cram_data->i = (size_t)-1; /* Current buffered trace */
        cram_data->ntr = 0; /* Number of buffered traces */
    } else {
        /* Local access mode via memory mapping */
        cram_data->mmaped = (char*)mmap (NULL, n + cram_data->offs,
					 PROT_READ, MAP_SHARED,
					 fileno (cram_data->stream), 0);
        if (cram_data->mmaped == MAP_FAILED)
            sf_error ("Data mmap failed: %s", strerror (errno));
        cram_data->data = (float*)(cram_data->mmaped + cram_data->offs);
    }
    if (cram_data->nc)
        sf_warning ("Using %d daemons for remote access to data", cram_data->nc);
    else
        sf_warning ("Using direct access to data file");

    return cram_data;
}

void sf_cram_data2_close (sf_cram_data2 cram_data)
/*< Destroy object >*/
{
    int i;
    size_t n = (size_t)cram_data->nt*(size_t)cram_data->n*(size_t)sizeof(float);

    if (cram_data->sockets) {
        for (i = 0; i < cram_data->nc; i++) {
            if (cram_data->sockets[i] != -1)
                close (cram_data->sockets[i]);
        }
        free (cram_data->sockets);
    }

    if (cram_data->nc) {
        sf_warning ("%g Gb of data was transfered through the network",
                    1e-9*(float)cram_data->nb);
        free (cram_data->trvals);
    } else
        munmap (cram_data->mmaped, n + cram_data->offs);
    free (cram_data);
}

bool sf_cram_data2_get_erefl (sf_cram_data2 cram_data)
/*< Return true if data is from exploding reflector modeling >*/
{
    return cram_data->erefl;
}

/* Return one sample from propely filtered trace according to requested slope */
static float sf_cram_data2_get_sample_from_trace (sf_cram_data2 cram_data, float *trace,
                                                  float trfact, float t, int *hits) {
    int it, k, ik1, ik2;
    float tf;

    tf = (t - cram_data->t0)/cram_data->dt;
    it = floorf (tf);
    tf -= (float)it;
    if (it < 0 || it >= (cram_data->nt - 1))
        return 0.0;

    if (cram_data->filter) {
        /* Filter length after Lumley-Claerbout-Bevc */
        k = 4.0*fabs (trfact) - 1.0;
        if (k < 1) k = 1;
        k = (k - 1)/2*2 + 1;
        ik1 = it - k/2 - 1;
        ik2 = it + k/2 + 1;
        if (ik1 < 0 || ik2 >= (cram_data->nt - 1))
            return 0.0;

        (*hits)++;
        return 1.0/((float)(k/2 + 1))*1.0/((float)(k/2 + 1))*(
                   2.0*(trace[it]*(1.0 - tf)
                        + trace[it + 1]*tf)
                  -1.0*(trace[ik1]*(1.0 - tf)
                        + trace[ik1 + 1]*tf)
                  -1.0*(trace[ik2]*(1.0 - tf)
                        + trace[ik2 + 1]*tf));
    } else {
        (*hits)++;
        return (trace[it]*(1.0 - tf)
                + trace[it + 1]*tf);
    }
}

/* In trace-by-trace mode, make sure that we have trace i currently buffered,
   return offset in the trace buffer */
static size_t sf_cram_data2_access_trace (sf_cram_data2 cram_data, size_t i) {
    size_t off;
    int len = 0, rc, is, elen;
    fd_set sset;
    struct timeval timeout;
    bool local = true;
    sf_cram_data_trreq trreq;

    /* Trace is not buffered */
#ifdef _OPENMP
#pragma omp critical
#endif
    if (i < cram_data->i || i >= (cram_data->i + cram_data->ntr)) {
        /* Data daemon index */
        is = (int)((float)i/cram_data->trd);
        if (cram_data->nc && cram_data->sockets[is] != -1) {
            len = 0; 
            /* Send trace request */
            trreq.id = rand ()*rand ();
            trreq.i = i;
            trreq.n = CRAM_TRBUF;
#ifdef TCP_CORK
            rc = 1;
            setsockopt (cram_data->sockets[is], SOL_TCP, TCP_CORK, &rc, sizeof(rc));
#endif
            while (len < sizeof(sf_cram_data_trreq)) {
                rc = send (cram_data->sockets[is], (const void*)(((unsigned char*)&trreq) + len),
                           sizeof(sf_cram_data_trreq) - len, MSG_NOSIGNAL);
                if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                    sf_warning ("Can not send trace request data for i=%lu, disconnecting from socket %d",
                                i, cram_data->sockets[is]);
                    close (cram_data->sockets[is]);
                    cram_data->sockets[is] = -1;
                    len = 0;
                    break;
                }
                if (rc > 0)
                    len += rc;
            }
#ifdef TCP_CORK
            rc = 0;
            setsockopt (cram_data->sockets[is], SOL_TCP, TCP_CORK, &rc, sizeof(rc));
#endif
            if (len) {
                FD_ZERO(&sset);
                FD_SET(cram_data->sockets[is], &sset);
                timeout.tv_sec = 3600;
                timeout.tv_usec = 0;
                len = 0;
                /* Wait for a response from the server */
                rc = select (cram_data->sockets[is] + 1, &sset, NULL, NULL, &timeout);
                if (0 == rc)
                    sf_warning ("Timeout reached for socket[%d]=%d while expecting a response for trace %lu",
                                is, cram_data->sockets[is], i);
                else if (rc < 0)
                    sf_error ("select() failed");
                else {
                    /* Receive reponse header from the server */
                    while (len < sizeof(sf_cram_data_trvals)) {
                        rc = recv (cram_data->sockets[is], (void*)(((unsigned char*)cram_data->trvals) + len),
                                   sizeof(sf_cram_data_trvals) - len, 0);
                        if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                            /* Connection terminated */
                            if (0 == rc)
                                sf_warning ("The server has closed connection for socket %d",
                                            cram_data->sockets[is]);
                            else
                                sf_warning ("Can not receive data for socket %d, disconnecting",
                                            cram_data->sockets[is]);
                            close (cram_data->sockets[is]);
                            len = 0;
                            break;
                        }
                        if (rc > 0)
                            len += rc;
                    }
                    /* Check correctness of the response header */
                    if (len) {
                        if (cram_data->trvals->id != trreq.id) {
                            sf_warning ("Received garbage from socket %d", cram_data->sockets[is]);
                            len = 0;
                        }
                        if (0 == cram_data->trvals->n) {
                            sf_warning ("Server replied that the trace index is out of bounds");
                            len = 0;
                        }
                    }
                    if (len) {
                        /* Expected length of traces */
                        elen = cram_data->trvals->n*cram_data->nt*sizeof(float);
                        if (cram_data->kmah)
                            elen *= 2;
                        /* Receive traces from the server */
                        len = 0;
                        while (len < elen) {
                            rc = recv (cram_data->sockets[is], (void*)(((unsigned char*)&cram_data->trvals->samples[0]) + len),
                                       elen - len, 0);
                            if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                                /* Connection terminated */
                                if (0 == rc)
                                    sf_warning ("The server has closed connection for socket %d",
                                                cram_data->sockets[is]);
                                else
                                    sf_warning ("Can not receive data for socket %d, disconnecting",
                                                cram_data->sockets[is]);
                                close (cram_data->sockets[is]);
                                len = 0;
                                break;
                            }
                            if (rc > 0)
                                len += rc;
                        }
                        if (len) {
                            local = false;
                            cram_data->ntr = cram_data->trvals->n;
                            cram_data->nb += (size_t)elen;
                        } else
                            sf_warning ("Only %d out of %d bytes received for trace %lu, expecting %n traces",
                                        len, elen, i, cram_data->trvals->n);
                    } /* Trace samples receive branch */
                }
            } /* Receive branch */
        } /* Send branch */
        if (local) {
            if (cram_data->nc && cram_data->nc && cram_data->sockets[is] != -1)
                sf_warning ("Failed to receive trace %lu from a remote server, using local access", i);
            if (i + (size_t)CRAM_TRBUF < cram_data->n)
                cram_data->ntr = CRAM_TRBUF;
            else
                cram_data->ntr = cram_data->n - i;
            off = (size_t)i*(size_t)cram_data->nt;
            if (cram_data->kmah)
                off *= (size_t)2;
            fseek (cram_data->stream, cram_data->offs + off*sizeof(float), SEEK_SET);
            fread (cram_data->data, sizeof(float), cram_data->kmah ? 2*cram_data->ntr*cram_data->nt
                                                                   : cram_data->ntr*cram_data->nt,
                   cram_data->stream);
        }
        cram_data->i = i;
    } 

    off = (i - cram_data->i)*(size_t)cram_data->nt;
    if (cram_data->kmah)
        off *= (size_t)2;
    return off;
}

float sf_cram_data2_get_sample (sf_cram_data2 cram_data, size_t i, float t,
                                float p, float dx, int kmah, int *hits)
/*< Get one sample for trace i, time t, slope p, and spatial sampling dx >*/
{
    size_t off = 0;
    float trf, w = 1.0;
    float *trace = NULL;

    if (cram_data->erefl) {
        t *= 0.5;
        kmah /= 2;
    }

    if (i >= cram_data->n)
        return 0.0;

    if (cram_data->nc) { /* Remote access mode */
        off = sf_cram_data2_access_trace (cram_data, i);
    } else { /* Local access mode with memory mapping */
        off = (size_t)i*(size_t)cram_data->nt;
        if (cram_data->kmah)
            off *= (size_t)2;
    }

    /* Choose largest trace factor for AA filter length */
    trf = fabs (p*dx/cram_data->dt);

    /* Choose trace according to KMAH index */
    if (cram_data->kmah) {
        kmah = kmah % 4; /* Enforce phase-shift periodicity */
    } else {
        kmah = 0;
    }

    switch (kmah) {
        case 0:
            trace = &cram_data->data[off];
            break;
        case 1: /* PI/2 */
            trace = &cram_data->data[off + (size_t)cram_data->nt];
            w = -1.0;
            break;
        case 2: /* PI */
            trace = &cram_data->data[off];
            w = -1.0;
            break;
        case 3: /* 3PI/4 */
            trace = &cram_data->data[off + (size_t)cram_data->nt];
            break;
    }

    return w*sf_cram_data2_get_sample_from_trace (cram_data, trace,
                                                  trf, t, hits);
}

