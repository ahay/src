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
#ifdef LINUX
#include <net/if.h>
#endif

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

#ifndef SOL_TCP
#define SOL_TCP IPPROTO_TCP
#endif

#include <rsf.h>

#include "einspline.h"
#include "esc_scgrid3.h"
#include "esc_helper.h"

/* This can be replaced with AF_INET_SDP(27) for Socket Direct Protocol */
#define ESC_SCD3_FAMILY AF_INET

/* Convert local domain name into an ASCII string with IP address */
#ifdef LINUX
static char *ip_linux = NULL;

static char* sf_escscd3_local_ip (int i) {
    int s;
    struct ifconf ifconf;
    struct ifreq ifr[50];
    struct sockaddr_in *s_in = NULL;
    int ifs;

    s = socket (AF_INET, SOCK_STREAM, 0);
    if (s < 0)
        sf_error ("socket() failed, errno=%d", errno);

    ifconf.ifc_buf = (char*)ifr;
    ifconf.ifc_len = sizeof(ifr);

    if (ioctl(s, SIOCGIFCONF, &ifconf) == -1)
        sf_error ("ioctl()[SIOCGIFCONF] failed, errno=%d", errno);
    ifs = ifconf.ifc_len/sizeof(ifr[0]);

    if (i >= ifs)
        sf_error ("Can not choose interface %d, only %d available", i, ifs);

    ip_linux = (char*)malloc(INET_ADDRSTRLEN);
    s_in = (struct sockaddr_in*)&ifr[i].ifr_addr;

    if (!inet_ntop (AF_INET, &s_in->sin_addr, ip_linux, INET_ADDRSTRLEN))
        sf_error ("inet_ntop() failed, errno=%d", errno);

    close (s);
    return ip_linux;
#else
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
#endif
}

/* Every servicing thread gets necessary data via this structure */
typedef struct {
    sf_esc_scgrid3_areq  *areqs;
    sf_esc_scgrid3_avals *avals;
    pthread_mutex_t      *smutex;
    int                  *njobs; /* Pointer to the job counter */
    size_t                ns, nr, ir, is; /* Send/receive statistics */
    int                   nd; /* Maximum number of requests to expect */
    int                   sd; /* Socket to read/write from/to */
    FILE                 *logf; /* Log file */
    int                   iab0, iab1, nab;
    multi_UBspline_3d_s  *scsplines;
} sf_escscd3_work;

static sf_escscd3_work* sf_escscd3_alloc_work (int nthreads, int ma, int mb) {
    sf_esc_scgrid3_areq *areqs;
    sf_esc_scgrid3_avals *avals;
    sf_escscd3_work* work;
    int i;

    work = (sf_escscd3_work*)sf_alloc (nthreads, sizeof(sf_escscd3_work));
    areqs = (sf_esc_scgrid3_areq*)sf_alloc (nthreads*ma*mb*SCGRID3_MAX_STENCIL,
                                            sizeof(sf_esc_scgrid3_areq));
    avals = (sf_esc_scgrid3_avals*)sf_alloc (nthreads*ma*mb*SCGRID3_MAX_STENCIL,
                                             sizeof(sf_esc_scgrid3_avals));
    for (i = 0; i < nthreads; i++) {
        work[i].areqs = &areqs[i*ma*mb*SCGRID3_MAX_STENCIL];
        work[i].avals = &avals[i*ma*mb*SCGRID3_MAX_STENCIL];
    }
    return work;
}

static void sf_escscd3_free_work (sf_escscd3_work* work) {
    free (work->avals);
    free (work->areqs);
    free (work);
}

/* This is a send/receive loop for each connection; performed in a separate thread */
static void* sf_escscd3_process_requests (void *ud) {
    sf_escscd3_work *data = (sf_escscd3_work*)ud;
    int iab;
    int i, len = 0, rc, n;
    fd_set fset;
    struct timeval timeout;

    /* Send/receive loop while connection exists */
    do {
        FD_ZERO(&fset);
        FD_SET(data->sd, &fset);
        timeout.tv_sec  = 86400;
        timeout.tv_usec = 0;
        rc = select (data->sd + 1, &fset, NULL, NULL, &timeout);
        if (rc <= 0) {
            if (rc == 0)
                fprintf (data->logf, "Timeout reached for socket %d, closing connection\n", data->sd);
            else
                fprintf (data->logf, "select() failed for socket %d, closing connection\n", data->sd);
            fflush (data->logf);
            close (data->sd);
            pthread_mutex_lock (data->smutex);
            (*data->njobs)--;
            pthread_mutex_unlock (data->smutex);
            return NULL;
        }
        /* Receive incoming data on this socket */
        len = 0;
        while (len < sizeof(sf_esc_scgrid3_areq)*data->nd) {
            /* Receive job request from the client */
            rc = recv (data->sd, (void*)(((unsigned char*)data->areqs) + len),
                       sizeof(sf_esc_scgrid3_areq)*data->nd - len, 0);
            if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                if (rc != 0)
                    fprintf (data->logf, "Connection was terminated for socket %d\n", data->sd);
                else
                    fprintf (data->logf, "Connection was closed for socket %d\n", data->sd);
                fprintf (data->logf, "%g GBytes received, %lu requests processed, average request %g KBytes for socket %d\n",
                         1e-9*(float)data->nr, data->ir, data->ir ? 1e-3*(float)data->nr/(float)data->ir : 0.0,
                         data->sd);
                fprintf (data->logf, "%g GBytes sent, %lu responses prepared, average response %g KBytes for socket %d\n",
                         1e-9*(float)data->ns, data->is, data->is ? 1e-3*(float)data->ns/(float)data->is : 0.0,
                         data->sd);
                fflush (data->logf);
                close (data->sd);
                pthread_mutex_lock (data->smutex);
                (*data->njobs)--;
                pthread_mutex_unlock (data->smutex);
                return NULL;
            }
            if (rc > 0)
                len += rc;
            if (0 == len % sizeof(sf_esc_scgrid3_areq)) /* A few complete requests have been received */
                break;
        }
        /* Process the received requests */
        n = len / sizeof(sf_esc_scgrid3_areq);
        data->ir++;
        data->nr += len;
        for (i = 0; i < n; i++) {
            data->avals[i].id = data->areqs[i].id;
            data->avals[i].ud1 = data->areqs[i].ud1;
            data->avals[i].ud2 = data->areqs[i].ud2;
            iab = data->areqs[i].iab;
            if (iab > data->iab1)
                iab -= data->nab;
            else if (iab < data->iab0)
                iab += data->nab;
            /* Interpolate */
            if (iab >= data->iab0 && iab <= data->iab1) {
                eval_multi_UBspline_3d_s (&data->scsplines[iab - data->iab0],
                                          data->areqs[i].y, data->areqs[i].x, data->areqs[i].z,
                                          data->avals[i].vals);
            } else {
                data->avals[i].vals[0] = SF_HUGE; /* Return error */
            }
        }
        /* Send the result back */
#ifdef TCP_CORK
        rc = 1;
        setsockopt (data->sd, SOL_TCP, TCP_CORK, &rc, sizeof(rc));
#endif
        len = 0;
        while (len < sizeof(sf_esc_scgrid3_avals)*n) {
            rc = send (data->sd, (const void*)(((unsigned char*)data->avals) + len),
                       sizeof(sf_esc_scgrid3_avals)*n - len, MSG_NOSIGNAL);
            if (rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) {
                fprintf (data->logf, "Connection was terminated for socket %d\n", data->sd);
                fflush (data->logf);
                close (data->sd);
                pthread_mutex_lock (data->smutex);
                (*data->njobs)--;
                pthread_mutex_unlock (data->smutex);
                return NULL;
            }
            if (rc > 0)
                len += rc;
        }
        data->is++;
        data->ns += len;
#ifdef TCP_CORK
        rc = 0;
        setsockopt (data->sd, SOL_TCP, TCP_CORK, &rc, sizeof(rc));
#endif
    } while (true);

    pthread_mutex_lock (data->smutex);
    (*data->njobs)--;
    pthread_mutex_unlock (data->smutex);
    return NULL;
}

#define MAX_BUF 1024

int main (int argc, char* argv[]) {
    size_t nc;
    off_t nsp;
    int ma, mb;
    int ith = 1, na, nb, nab, ncv, icpu, ncpu, tout, tdel;
    int i, iab, iab0, iab1, port, nthreads, tmpfile = 0;
#ifdef LINUX
    int inet = 0;
#endif
    multi_UBspline_3d_s *scsplines = NULL;
    sf_file in, scgrid = NULL, out;
    FILE *logf;
    char sbuffer[MAX_BUF];
    char *str = NULL;
    /* Fork variables */
    pid_t pid, sid;
    /* Server network variables */
    char *ip = NULL;
    int rc, on = 1, ijob, bsiz;
    int listen_sd=0, new_sd, njobs = 0;
    struct sockaddr_in serv_addr, client_addr;
    struct timeval timeout;
    fd_set sset;
    socklen_t clen;
    /* Threading */
    sf_escscd3_work *rwork, *qjobs, *old_qjobs;
    pthread_t pthread;
    pthread_mutex_t smutex;

    sf_init (argc, argv);

    in = sf_input ("in");
    /* Angular grid map */

    out = sf_output ("out");
    /* Dummy output */

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

#ifdef LINUX
    if (!sf_getint ("inet", &inet)) inet = 1;
    /* Network interface index */
#endif

    if (ith) {
#ifdef LINUX
        if ((ip = sf_escscd3_local_ip (inet)))
#else
        if ((ip = sf_escscd3_local_ip ()))
#endif
        {
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
    /* How many azimuth angles to expect per request */
    if (!sf_getint ("mb", &mb)) mb = 20;
    /* How many inclination angles to expect per request */
    if (!sf_getint ("nthreads", &nthreads)) nthreads = 2*ncpu;
    /* Number of threads (connections) per daemon */
    if (!sf_getint ("timeout", &tout)) tout = 10;
    /* Inactivity time before shutdown (mins) */
    if (!sf_getint ("tdelay", &tdel)) tdel = 0;
    /* Time delay before accessing data, tdel*icpu (secs) */

    if (ith && 0 == (icpu % ith))
        sf_warning ("Running no more than %d threads on %s [CPU %d]", nthreads, ip, icpu);

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

    if (nab >= na*nb)
        nab = na*nb;
    /* Number of daemons for a full coverage */
    ncv = na*nb/nab;
    if ((na*nb) % nab)
        sf_error ("na*nb should be divisible by nab");
    if (ith && (ncpu % ncv))
        sf_error ("Non-integer number of coverages");
    if (ncpu*nab < na*nb)
        sf_error ("Incomplete angle coverage");
    sf_putint (out, "Nab", nab);
    sf_putint (out, "Ncv", ncv);

    /* Determine how much data from the supercell grid to hold per CPU */
    if (ith && 0 == (icpu % ith)) {
        iab0 = ((icpu / ith)*nab) % (na*nb);
        iab1 = iab0 + nab - 1;
        sf_warning ("Serving angular patch from iab=%d to iab=%d [CPU %d]",
                    iab0, iab1, icpu);
    } else {
        /* No daemon for this CPU */
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
            sf_error ("fork() failed, errno=%d", errno);
        if (0 == pid) {
            if (lockf (tmpfile, F_LOCK, 0) == -1)
                abort();
        }
        else
            sleep (1);
    } else
        pid = 1;

    /* Read spline coefficients for the determined range of angles */
    if (ith && 0 == (icpu % ith) && 0 == pid) {
        /* Stream socket for incoming connections on */
        listen_sd = socket (ESC_SCD3_FAMILY, SOCK_STREAM, 0);
        if (listen_sd < 0)
            sf_error ("socket() failed [CPU %d], errno=%d", icpu, errno);
/*      
        new_sd = connect (listen_sd, (struct sockaddr *)&serv_addr, sizeof(serv_addr));
        if (0 == new_sd) {
            sf_warning ("Daemon is already running [CPU %d]", icpu);
            close (new_sd);
            close (listen_sd);
            lockf (tmpfile, F_ULOCK, 1);
            close (tmpfile);
            free (str);
            return 0;
        }
*/
        sleep ((tdel + 1)*(icpu/ith));
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
        sf_warning ("Exiting parent process [CPU %d]", icpu);
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
        sf_error ("setsid() failed [CPU %d], errno=%d", icpu, errno);
    /* Change to the root directory to prevent locking the current one */
/*
    if ((chdir ("/")) < 0)
        sf_error ("chdir() failed [CPU %d]", icpu);
*/

    /*************************************/
    /* Server part, use non-blocking I/O */
    /*************************************/

    /* Allow socket descriptor to be reuseable */
    if (setsockopt (listen_sd, SOL_SOCKET, SO_REUSEADDR,
                    (char *)&on, sizeof(on)) < 0) {
        close (listen_sd);
        sf_error ("setsockopt() failed [CPU %d], errno=%d", icpu, errno);
    }
    /* Set socket to be non-blocking; all of the sockets for
       the incoming connections will also be non-blocking */
    if (ioctl (listen_sd, FIONBIO, (char *)&on) < 0) {
        close (listen_sd);
        sf_error ("ioctl() failed [CPU %d], errno=%d", icpu, errno);
    }

    /* Bind the socket */
    memset (&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port = htons(port);
    if (bind (listen_sd,
              (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        close (listen_sd);
        sf_error ("bind() failed [CPU %d], errno=%d", icpu, errno);
    }

    /* Set the listen back log */
    if (listen (listen_sd, nthreads) < 0) {
        close (listen_sd);
        sf_error ("listen() failed [CPU %d], errno=%d", icpu, errno);
    }

    pthread_mutex_init (&smutex, NULL);

    /* Open log file */
    snprintf (sbuffer, MAX_BUF, "sfescscd3_%s_%d.log", ip, icpu/ith);
    logf = fopen (sbuffer, "w+");
    fprintf (logf, "Listening on %s [CPU %d]\n", ip, icpu);
    fprintf (logf, "Servicing angle patch %d - %d\n", iab0, iab1);
    fflush (logf);
    sf_warning ("Log file is %s [CPU %d]", sbuffer, icpu);

    /* Buffer for job requests */
    old_qjobs = NULL;
    qjobs = sf_escscd3_alloc_work (nthreads, ma, mb);
    /* Position in the buffer for the next client */
    ijob = 0;
    clen = sizeof(client_addr);

    /* Release the child before starting the server loop */
    if (lockf (tmpfile, F_ULOCK, 1) == -1)
        abort();
    close(tmpfile);
    free (str);
    fclose (stderr);

    /* Wait for incoming connections */
    do {
        FD_ZERO(&sset);
        FD_SET(listen_sd, &sset);
        timeout.tv_sec  = tout*60;
        timeout.tv_usec = 0;
        /* Wait for incoming connections */
        rc = select (listen_sd + 1, &sset, NULL, NULL, &timeout);
        if (rc < 0) {
            fprintf (logf, "select() failed\n");
            fflush (logf);
            return -1;
        }

        /* Check to see if the timeout has expired */
        if (0 == rc) {
            pthread_mutex_lock (&smutex);
            rc = njobs;
            pthread_mutex_unlock (&smutex);            
            if (0 == rc) {
                tout = -1;
                break;
            }
            fprintf (logf, "select() timeout, work is still in progress, postponing shutdown\n");
            fflush (logf);
            continue;
        }

        /* Accept the incoming connection */
        new_sd = accept (listen_sd, (struct sockaddr*)&client_addr, &clen);
        if (new_sd < 0) {
            if (errno != EWOULDBLOCK && errno != EAGAIN) {
                fprintf (logf, "accept() failed\n");
                fflush (logf);
                return -1;
            }
            break;
        }
        ip = inet_ntoa (client_addr.sin_addr);
        on = 1;
        if (ioctl (new_sd, FIONBIO, (char *)&on) < 0) {
            fprintf (logf, "ioctl() failed for a new socket\n");
            fflush (logf);
            close (new_sd);
            break;
        }
        on = 1;
        if (setsockopt (new_sd, SOL_SOCKET, SO_REUSEADDR,
                        (char *)&on, sizeof(on)) < 0) {
            fprintf (logf, "setsockopt() failed for a new socket\n");
            fflush (logf);
            close (new_sd);
            break;
        }
        on = 1;
#ifndef TCP_CORK
        if (setsockopt (new_sd, SOL_TCP, TCP_NODELAY, &on, sizeof(on)) < 0) { 
            fprintf (logf, "Can not set TCP_NODELAY for a new connection\n");
            close (new_sd);
            break;
        }
#endif
#ifdef SO_NOSIGPIPE
        on = 1;
        if (setsockopt (new_sd, SOL_SOCKET, SO_NOSIGPIPE, &on, sizeof(on)) < 0) {
            fprintf (logf, "Can not set SO_NOSIGPIPE for a new connection\n");
            fprintf (logf, "Rejecting connection from %s\n", ip);
            fflush (logf);
            close (new_sd);
            break;
        }
#endif
        bsiz = sizeof(sf_esc_scgrid3_areq)*ma*mb*SCGRID3_MAX_STENCIL;
        if (setsockopt (new_sd, SOL_SOCKET, SO_RCVBUF, &bsiz, sizeof(int)) < 0) {
            fprintf (logf, "Can not set SO_RCVBUF for a new connection\n");
            close (new_sd);
            break;
        }
        bsiz = sizeof(sf_esc_scgrid3_avals)*ma*mb*SCGRID3_MAX_STENCIL;
        if (setsockopt (new_sd, SOL_SOCKET, SO_SNDBUF, &bsiz, sizeof(int)) < 0) {
            fprintf (logf, "Can not set SO_SNDBUF for a new connection\n");
            close (new_sd);
            break;
        }

        /* Connection with the new client has been established,
           prepare data for processing its requests */
        rwork = &qjobs[ijob];                    
        rwork->sd = new_sd;
        rwork->iab0 = iab0;
        rwork->iab1 = iab1;
        rwork->nab = na*nb;
        rwork->scsplines = scsplines;
        rwork->logf = logf;
        rwork->nd = ma*mb*SCGRID3_MAX_STENCIL;
        rwork->smutex = &smutex;
        rwork->njobs = &njobs;
        rwork->is = 0;
        rwork->ir = 0;
        rwork->ns = 0;
        rwork->nr = 0;

        pthread_mutex_lock (&smutex);
        if (njobs < nthreads) {
            if (pthread_create (&pthread, NULL,
                                sf_escscd3_process_requests, (void*)rwork) != 0) {
                fprintf (logf, "pthread_create() failed\n");
                close (new_sd);
                pthread_mutex_unlock (&smutex);
                break;
            }
            njobs++;
            ijob++;
            /* Find a free block for the new client incoming requests */
            if (ijob == nthreads) {
                if (old_qjobs)
                    sf_escscd3_free_work (old_qjobs);
                old_qjobs = qjobs;
                qjobs = sf_escscd3_alloc_work (nthreads, ma, mb);
                ijob = 0;
            }
            fprintf (logf, "Accepted client from %s, socket %d\n", ip, new_sd);
            fflush (logf);
        } else {
            fprintf (logf, "Reached maximum allowed number of threads, rejecting connection from %s\n", ip);
            close (new_sd);
        }
        pthread_mutex_unlock (&smutex);
    } while (true);

    if (tout < 0)
        fprintf (logf, "Shutting down on timeout\n");
    fflush (logf);
    close (listen_sd);

    /* Clean up */
    pthread_mutex_destroy (&smutex);
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

