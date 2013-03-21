/* Daemon for distributed computation of stitched escape solutions in supercells in 3-D. */
/*
  Copyright (C) 2013 University of Texas at Austin

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

#ifdef PTHREADS

#include <pthread.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/ioctl.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>

#include <rsf.h>

#include "einspline.h"
#include "esc_scgrid3.h"
#include "esc_helper.h"

/* Convert local domain name into an ASCII string with IP address */
static char* sf_escscd3_local_ip () {
    char hostname[1024];
    struct hostent *he;
    struct in_addr **addr_list;

    if (gethostname (hostname, 1024) != 0 ||
        (he = gethostbyname (hostname)) == NULL)
        return NULL;

    addr_list = (struct in_addr **)he->h_addr_list;
    if (addr_list[0] != NULL)
        return (inet_ntoa (*addr_list[0]));

    return NULL;
}

/* Chunk of work to be added to the queue */
typedef struct {
    sf_esc_scgrid3_areq  *areqs;
    sf_esc_scgrid3_avals *avals;
    int                   n; /* Number of requests */
    int                   sd; /* Socket to write to */
    pthread_mutex_t      *smutex; /* Mutex for serial access to the socket */
    int                   iab0, iab1, nab;
    multi_UBspline_3d_s  *scsplines;
} sf_escscd3_work;

static sf_escscd3_work* sf_escscd3_alloc_work (int mjobs, int ma, int mb) {
    sf_esc_scgrid3_areq *areqs;
    sf_esc_scgrid3_avals *avals;
    sf_escscd3_work* work;
    int i;

    work = (sf_escscd3_work*)sf_alloc (mjobs, sizeof(sf_escscd3_work));
    areqs = (sf_esc_scgrid3_areq*)sf_alloc (mjobs*ma*mb, sizeof(sf_esc_scgrid3_areq));
    avals = (sf_esc_scgrid3_avals*)sf_alloc (mjobs*ma*mb, sizeof(sf_esc_scgrid3_avals));
    for (i = 0; i < mjobs; i++) {
        work[i].areqs = &areqs[i*ma*mb];
        work[i].avals = &avals[i*ma*mb];
    }
    return work;
}

static void sf_escscd3_free_work (sf_escscd3_work* work) {
    free (work->avals);
    free (work->areqs);
    free (work);
}

/* Elementary task to be peformed by the workers */
static void sf_escscd3_extract_point (void *ud) {
    sf_escscd3_work *data = (sf_escscd3_work*)ud;
    int iab;
    int i, len = 0, rc;

    for (i = 0; i < data->n; i++) {
        data->avals[i].id = data->areqs[i].id;
        data->avals[i].ud1 = data->areqs[i].ud1;
        data->avals[i].ud2 = data->areqs[i].ud2;
        iab = data->areqs[i].iab;
        if (iab > data->iab1)
            iab -= data->nab;
        else if (iab < data->iab0)
            iab += data->nab;
        /* Interpolate */
        if (iab >= data->iab0 && iab <= data->iab1)
            eval_multi_UBspline_3d_s (&data->scsplines[iab - data->iab0],
                                      data->areqs[i].y, data->areqs[i].x, data->areqs[i].z,
                                      data->avals[i].vals);
        else
            data->avals[i].vals[0] = SF_HUGE; /* Return error */
    }
    /* Send the result back */
    pthread_mutex_lock (data->smutex);
    while (len < sizeof(sf_esc_scgrid3_avals)*data->n) {
        rc = send (data->sd, (const void*)(((unsigned char*)data->avals) + len),
                   sizeof(sf_esc_scgrid3_avals)*data->n - len, 0);
        if (rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK)
            break; /* The connection is gone */
        if (rc > 0)
            len += rc;
        else
            break;
    }
    pthread_mutex_unlock (data->smutex);
}

#define MAX_BUF 1024

int main (int argc, char* argv[]) {
    size_t nc;
    off_t nsp;
    int ma, mb;
    int ith = 1, n1, na, nb, nab, icpu, ncpu, bcpu, wcpu, tout;
    int i, iab, iab0, iab1, port, nthreads, tmpfile = 0;
    multi_UBspline_3d_s *scsplines = NULL;
    sf_file in, scgrid = NULL, out;
    FILE *logf;
    char sbuffer[MAX_BUF];
    char *str = NULL;
    /* Fork variables */
    pid_t pid, sid;
    /* Server network variables */
    char *ip = NULL;
    int len, rc, on = 1, ijob, mjobs;
    int listen_sd, max_sd, new_sd, desc_ready;
    bool close_conn = false;
    struct sockaddr_in serv_addr, client_addr;
    struct timeval timeout;
    fd_set master_set, working_set;
    socklen_t clen;
    /* Thread pool */
    sf_thpool tpool;
    sf_escscd3_work *rwork, *qjobs, *old_qjobs;
    pthread_mutex_t smutex[FD_SETSIZE];

    sf_init (argc, argv);

    in = sf_input ("in");
    /* Angular grid map */

    out = sf_output ("out");
    /* Dummy output */

    if (!sf_histint (in, "n1", &n1)) sf_error ("Need n1=");

    if (!sf_histint (in, "icpu", &icpu)) icpu = 0;
    /* Current CPU number */
    if (!sf_histint (in, "ncpu", &ncpu)) ncpu = 1;
    /* Total number of CPUs */

    if (!sf_getint ("nab", &nab)) nab = 1;
    /* Number of angular blocks to keep in memory per daemon */
    if (!sf_getint ("port", &port)) port = 29542;
    /* TCP port for listening */
    if (!sf_getint ("ith", &ith)) ith = 0;
    /* Make every ith process a daemon */

    memset (&serv_addr, 0, sizeof (serv_addr));

    if (ith) {
        if ((ip = sf_escscd3_local_ip ())) {
            if (0 == icpu % ith)
                sf_warning ("Assuming IP address %s [CPU %d]", ip, icpu);
            serv_addr.sin_family = AF_INET; /* Internet address family */
            serv_addr.sin_addr.s_addr = inet_addr (ip);   /* Server IP address */
            serv_addr.sin_port = htons ((uint16_t)port); /* Server port */
        } else
            sf_error ("Can not determine local IP address, shutting down [CPU %d]", icpu);
    } else
        sf_warning ("No processes are selected to run as daemons, shutting down [CPU %d]", icpu);

    if (ith && icpu % ith)
        sf_warning ("Making room for the daemon on CPU %d, shutting down [CPU %d]", (icpu/ith)*ith, icpu);

    if (!sf_getint ("ma", &ma)) ma = 20;
    /* How many azimuth angles to expect per job */
    if (!sf_getint ("mb", &mb)) mb = 20;
    /* How many inclination angles to expect per job */
    if (!sf_getint ("mjobs", &mjobs)) mjobs = ncpu*10;
    /* Maximum number of jobs to hold in the queue */
    if (!sf_getint ("nthreads", &nthreads)) nthreads = ith;
    /* Number of threads per daemon */
    if (!sf_getint ("timeout", &tout)) tout = 10;
    /* Inactivity time before shutdown (mins) */

    if (ith && 0 == (icpu % ith))
        sf_warning ("Running %d threads on %s [CPU %d]", nthreads, ip, icpu);

    if (!sf_getstring ("scgrid")) sf_error ("Need scgrid=");
    /* Grid of supercells of local escape solutions */
    scgrid = sf_input ("scgrid");

    if (!sf_histlargeint (scgrid, "n1", &nsp)) sf_error ("No n1= in supercell file");
    /* Size of one angular direction */
    if (!sf_histint (scgrid, "Nb", &nb)) sf_error ("No Nb= in supercell file");
    /* Number of inclination angles */
    if (!sf_histint (scgrid, "Na", &na)) sf_error ("No Na= in supercell file");
    /* Number of azimuth angles */

    sf_settype (out, SF_UCHAR);
    sf_putlargeint (out, "n1", sizeof (serv_addr) + 2*sizeof (int));
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Supercell grid distributed access info");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", 1);
    sf_putfloat (out, "o2", 0.0);
    sf_putfloat (out, "d2", 1.0);
    sf_putstring (out, "label2", "");
    sf_putstring (out, "unit2", "");
    sf_putint (out, "n3", 1);
    sf_putfloat (out, "o3", 0.0);
    sf_putfloat (out, "d3", 1.0);
    sf_putstring (out, "label3", "");
    sf_putstring (out, "unit3", "");

    sf_putint (out, "Ndaemon", ith ? ncpu/ith : 0);
    sf_putint (out, "Port", port);
    sf_putint (out, "Ma", ma);
    sf_putint (out, "Mb", mb);
    sf_putstring (out, "Remote", ith ? "y" : "n");

    if (ith && nab*(ncpu/ith) < na*nb)
        sf_error ("Incomplete angle coverage; increase nab= or number of CPUs");
    /* Determine how much data from the supercell grid to hold per CPU */
    if (ith && 0 == (icpu % ith)) {
        if (nab >= na*nb) {
            iab0 = 0;
            iab1 = na*nb - 1;
        } else {
            /* Find out displacement */
            wcpu = (int)((float)(na*nb/(ncpu/ith)) + 1.0);
            bcpu = na*nb - (ncpu/ith)*(wcpu - 1);
            if ((icpu/ith) < bcpu)
                iab0 = (icpu/ith)*wcpu;
            else
                iab0 = bcpu*wcpu + (icpu/ith - bcpu)*(wcpu - 1);
            iab1 = iab0 + nab/2;
            if (nab % 2)
                iab0 -= nab/2;
            else
                iab0 -= (nab/2 - 1);
        }
        sf_warning ("Serving angular patch from iab=%d to iab=%d [CPU %d]",
                    iab0, iab1, icpu);
    } else {
        iab0 = nab;
        iab1 = -nab;
    }

    sf_ucharwrite ((unsigned char*)&serv_addr, sizeof (serv_addr), out);
    sf_ucharwrite ((unsigned char*)&iab0, sizeof (int), out);
    sf_ucharwrite ((unsigned char*)&iab1, sizeof (int), out);

    if (ith && 0 == (icpu % ith)) {
        /* Temporary file for synchronizing forked processes */
        str = strdup ("/tmp/sfescscd3.XXXXXX");
        tmpfile = mkstemp (str);
        /* Daemonize */
        pid = fork ();
        if (pid < 0)
            sf_error ("fork() failed");
        if (0 == pid)
            lockf (tmpfile, F_LOCK, 0);
        else
            sleep (1);
    } else
        pid = 1;

    /* Read spline coefficients for the determined range of angles */
    if (ith && 0 == (icpu % ith) && 0 == pid) {
        nc = 0;
        scsplines = (multi_UBspline_3d_s*)sf_alloc ((size_t)(iab1 - iab0 + 1),
                                                    sizeof(multi_UBspline_3d_s));
        for (iab = iab0; iab <= iab1; iab++) {
            i = iab;
            if (i < 0)
                i += na*nb;
            else if (i >= na*nb)
                i -= na*nb;
            sf_seek (scgrid, (off_t)i*(off_t)nsp, SEEK_SET);
            sf_ucharread ((unsigned char*)&scsplines[iab - iab0], sizeof(multi_UBspline_3d_s), scgrid);
            scsplines[iab - iab0].coefs = sf_floatalloc (scsplines[iab - iab0].nc/sizeof(float));
            sf_ucharread ((unsigned char*)scsplines[iab - iab0].coefs, scsplines[iab - iab0].nc, scgrid);
            nc += scsplines[iab - iab0].nc;
        }
        sf_warning ("%g Mb of spline coefficients buffered [CPU %d]", 1e-6*(float)nc, icpu);
        sf_warning ("Running as a daemon now [CPU %d]", icpu);
    } 

    if (tmpfile && pid != 0) {
        /* Wait for the parent to complete data preparation */
        while (lockf (tmpfile, F_LOCK, 1))
            sleep (1);
        sf_warning ("Exiting parent process");
        close (tmpfile);
        unlink (str);
        free (str);
    }

    sf_fileclose (scgrid);
    sf_fileclose (in);
    sf_fileclose (out);

    if (pid > 0)
        return 0;

    /* Change the file mode mask */
    umask (0);
    /* Create a new SID for the child process */
    sid = setsid ();
    if (sid < 0)
        sf_error ("setsid() failed");
    /* Change to the root directory to prevent locking the current one */
/*
    if ((chdir ("/")) < 0)
        sf_error ("chdir() failed");
*/
    /* Initialize workers */
    tpool = sf_thpool_init (nthreads, 2*mjobs);

    /*************************************/
    /* Server part, use non-blocking I/O */
    /*************************************/

    /* Stream socket for incoming connections on */
    listen_sd = socket (AF_INET, SOCK_STREAM, 0);
    if (listen_sd < 0)
        sf_error ("socket() failed");

    /* Allow socket descriptor to be reuseable */
    if (setsockopt (listen_sd, SOL_SOCKET, SO_REUSEADDR,
                    (char *)&on, sizeof(on)) < 0) {
        close (listen_sd);
        sf_error ("setsockopt() failed");
    }

    /* Set socket to be non-blocking; all of the sockets for
       the incoming connections will also be non-blocking */
    if (ioctl (listen_sd, FIONBIO, (char *)&on) < 0) {
        close (listen_sd);
        sf_error ("ioctl() failed");
    }

    /* Bind the socket */
    memset (&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family      = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port        = htons(port);
    if (bind (listen_sd,
              (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        close (listen_sd);
        sf_error ("bind() failed");
    }

    /* Set the listen back log */
    if (listen (listen_sd, 32) < 0) {
        close (listen_sd);
        sf_error ("listen() failed");
    }

    /* Initialize the master set of file descriptors */
    FD_ZERO(&master_set);
    max_sd = listen_sd;
    FD_SET(listen_sd, &master_set);

    /* Open log file */
    snprintf (sbuffer, MAX_BUF, "sfescscd3_%s_%d.log", ip, icpu/ith);
    logf = fopen (sbuffer, "w+");
    fprintf (logf, "Listening on %s\n", ip);
    fprintf (logf, "Servicing angle patch %d - %d\n", iab0, iab1);
    fflush (logf);
    sf_warning ("Log file is %s [CPU %d]", sbuffer, icpu);

    /* Release the child before starting the server loop */
    lockf (tmpfile, F_ULOCK, 1);
    close (tmpfile);
    free (str);
    fclose (stderr);

    /* Buffer for job requests */
    old_qjobs = NULL;
    qjobs = sf_escscd3_alloc_work (mjobs, ma, mb);
    /* Position in the buffer for the next job */
    ijob = 0;
    clen = sizeof(client_addr);

    /* Loop waiting for incoming connects or for incoming data */
    do {
        /* Die after a certain period of inactivity */
        timeout.tv_sec  = tout*60;
        timeout.tv_usec = 0;
        memcpy (&working_set, &master_set, sizeof(master_set));
        /* Poll sockets for incoming data */
        rc = select (max_sd + 1, &working_set, NULL, NULL, &timeout);
        if (rc < 0) {
            fprintf (logf, "select() failed");
            fflush (logf);
            return -1;
        }

        /* Check to see if the 5 minute time out expired */
        if (rc == 0)
            break;

        /* Check to see if this is the listening socket */
        if (rc && FD_ISSET (listen_sd, &working_set)) {
            /* Accept the incoming connection */
            new_sd = accept (listen_sd, (struct sockaddr*)&client_addr, &clen);
            if (new_sd < 0) {
                if (errno != EWOULDBLOCK && errno != EAGAIN) {
                    fprintf (logf, "accept() failed");
                    fflush (logf);
                    return -1;
                }
                break;
            }
            ip = inet_ntoa (client_addr.sin_addr);
            if (new_sd >= FD_SETSIZE) {
                fprintf (logf, "Maximum number of connections (FD_SETSIZE) is exceeded\n");
                fprintf (logf, "Rejecting connection from %s\n", ip);
                fflush (logf);
                close (new_sd);
                break;
            }
            on = 1;
            if (ioctl (new_sd, FIONBIO, (char *)&on) < 0) {
                fprintf (logf, "ioctl() failed for a new socket\n");
                fflush (logf);
                close (new_sd);
                break;
            }
            on = 1;
            if (setsockopt (new_sd, SOL_TCP, TCP_NODELAY, &on, sizeof(on)) < 0) {
                fprintf (logf, "Can not set TCP_NODELAY for a new connection\n");
                fprintf (logf, "Rejecting connection from %s\n", ip);
                fflush (logf);
                close (new_sd);
                break;
            }
            fprintf (logf, "Accepted client from %s, socket %d\n", ip, new_sd);
            fflush (logf);
            /* Add the new incoming connection to the master read set */
            FD_SET(new_sd, &master_set);
            if (new_sd > max_sd)
                max_sd = new_sd;
            pthread_mutex_init (&smutex[new_sd], NULL);
            FD_CLR(listen_sd, &working_set);
            rc--;
        }

        /* One or more descriptors are readable */
        desc_ready = rc;
        for (i = 0; i <= max_sd && desc_ready; i++) {
            /* Check to see if this descriptor is ready */
            if (!FD_ISSET (i, &working_set))
                continue;
            desc_ready--;
            rwork = &qjobs[ijob];                    
            rwork->sd = i;
            rwork->iab0 = iab0;
            rwork->iab1 = iab1;
            rwork->nab = na*nb;
            rwork->scsplines = scsplines;
            rwork->smutex = &smutex[i];
            close_conn = false;
            /* Receive incoming data on this socket */
            len = 0;
            while (len < sizeof(sf_esc_scgrid3_areq)*ma*mb) {
                /* Receive job request from the client */
                rc = recv (i, (void*)(((unsigned char*)rwork->areqs) + len),
                           sizeof(sf_esc_scgrid3_areq)*ma*mb - len, 0);
                if (rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) {
                    /* Connection terminated */
                    fprintf (logf, "Connection with socket %d is terminated\n", i);
                    fflush (logf);
                    close_conn = true;
                    break;
                }
                if (0 == rc) {
                    /* Normal shutdown */
                    fprintf (logf, "Socket %d has been closed by the remote party\n", i);
                    fflush (logf);
                    close_conn = true;
                    break;
                }
                if (rc > 0)
                    len += rc;
                if ((rc < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) || /* Nothing to receive */
                    0 == len % sizeof(sf_esc_scgrid3_areq)) /* A few complete requests have been received */
                    break;
            }
            if (close_conn) {
                /* Clean up after disconnect */
                pthread_mutex_lock (&smutex[i]);
                close (i);
                pthread_mutex_unlock (&smutex[i]);
                FD_CLR(i, &master_set);
                if (i == max_sd) {
                    while (FD_ISSET(max_sd, &master_set) == false)
                        max_sd -= 1;
                }
                pthread_mutex_destroy (&smutex[i]);
                len = 0;
            }
            if (len > 0) { /* Check if received a complete set of requests */
                if (len % sizeof(sf_esc_scgrid3_areq)) {
                    fprintf (logf, "Partial receive from socket %d\n", i);
                    fflush (logf);
                } else { /* Otherwise - add to the work queue */
                    rwork->n = len / sizeof(sf_esc_scgrid3_areq);
                    if (SF_THPOOL_QUEUE_FULL == /* Add requests to the work queue */
                        sf_thpool_add_work (tpool, sf_escscd3_extract_point, (void*)rwork)) {
                        fprintf (logf, "The queue is full, dropping request from socket %d\n", i);
                    } else {
                        ijob++;
                        /* Find a free block for the next job request */
                        if (ijob == mjobs) {
                            if (old_qjobs && sf_thpool_jobsn (tpool) <= mjobs)
                                sf_escscd3_free_work (old_qjobs);
                            old_qjobs = qjobs;
                            qjobs = sf_escscd3_alloc_work (mjobs, ma, mb);
                            ijob = 0;
                        }
                    }
                }
            }
        } /* Loop through selectable descriptors */
    } while (true);

    fprintf (logf, "%d jobs are left in the queue\n", sf_thpool_jobsn (tpool));
    fprintf (logf, "Shutting down on timeout\n");
    fflush (logf);

    /* Shutdown all connections */
    for (i=0; i <= max_sd; i++) {
        if (FD_ISSET(i, &master_set)) {
            pthread_mutex_destroy (&smutex[i]);
            close (i);
        }
    }
    /* Clean up */
    sf_thpool_close (tpool);
    if (old_qjobs)
        sf_escscd3_free_work (old_qjobs);
    sf_escscd3_free_work (qjobs);

    for (iab = iab0; iab <= iab1; iab++) {
        free (scsplines[iab - iab0].coefs);
    }
    free (scsplines);
    fclose (logf);

    return 0;
}

#else /* PTHREADS */

#include <rsf.h>

int main (int argc, char* argv[]) {
    fprintf (stderr, "sfescscd3 can not work without pthreads; reconfigure Madagascar with pthreads support and recompile\n");
    return -1;
}
#endif

