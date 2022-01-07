/* Daemon for distributed storage of prestack data for angle migration. */
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

#include "cram_data2.h"

/* This can be replaced with AF_INET_SDP(27) for Socket Direct Protocol */
#define CRAM_DATA_FAMILY AF_INET

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

/* Every connection servicing thread gets necessary data via this structure */
typedef struct {
    sf_cram_data_trreq   *trreq;
    sf_cram_data_trvals  *trvals;
    pthread_mutex_t      *smutex;
    int                  *njobs; /* Pointer to the job counter */
    int                   sd; /* Socket to read/write from/to */
    FILE                 *logf; /* Log file */
    size_t                i0, i1; /* Min/max trace indices */
    int                   nt; /* Number of time samples in each trace */
    bool                  kmah;
    float                *traces;
} sf_cramdd_work;

static sf_cramdd_work* sf_cramdd_alloc_work (int nthreads) {
    sf_cram_data_trreq *trreqs;
    sf_cram_data_trvals *trvals;
    sf_cramdd_work* work;
    int i;

    work = (sf_cramdd_work*)sf_alloc (nthreads, sizeof(sf_cramdd_work));
    trreqs = (sf_cram_data_trreq*)sf_alloc (nthreads, sizeof(sf_cram_data_trreq));
    trvals = (sf_cram_data_trvals*)sf_alloc (nthreads, sizeof(sf_cram_data_trvals));
    for (i = 0; i < nthreads; i++) {
        work[i].trreq = &trreqs[i];
        work[i].trvals = &trvals[i];
    }
    return work;
}

static void sf_cramdd_free_work (sf_cramdd_work* work) {
    free (work->trvals);
    free (work->trreq);
    free (work);
}

/* This is a send/receive loop for each connection; performed in a separate thread */
static void* sf_cramdd_process_requests (void *ud) {
    sf_cramdd_work *data = (sf_cramdd_work*)ud;
    int len = 0, rc, elen;
    float *traces = NULL;
    fd_set fset;
    struct timeval timeout;

    /* Send/receive loop while connection exists */
    do {
        FD_ZERO(&fset);
        FD_SET(data->sd, &fset);
        timeout.tv_sec = 86400;
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
        while (len < sizeof(sf_cram_data_trreq)) {
            /* Receive job request from the client */
            rc = recv (data->sd, (void*)(((unsigned char*)data->trreq) + len),
                       sizeof(sf_cram_data_trreq) - len, 0);
            if ((rc < 0 && errno != EAGAIN && errno != EWOULDBLOCK) || 0 == rc) {
                if (rc != 0)
                    fprintf (data->logf, "Connection was terminated for socket %d\n", data->sd);
                else
                    fprintf (data->logf, "Connection was closed for socket %d\n", data->sd);
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
        /* Process the received request */
        data->trvals->id = data->trreq->id;
        if (data->trreq->i < data->i0 || data->trreq->i > data->i1)
            /* Request is out of range */
            data->trvals->n = 0;
        else if (data->trreq->i + data->trreq->n > (data->i1 + (size_t)1))
            /* Return fewer traces than requested */
            data->trvals->n = data->i1 - data->trreq->i + (size_t)1;
        else
            data->trvals->n = data->trreq->n;
        if (data->trvals->n) {
            traces = data->kmah ? &data->traces[(data->trreq->i - data->i0)*(size_t)data->nt*(size_t)2]
                                : &data->traces[(data->trreq->i - data->i0)*(size_t)data->nt];
        }
        /* Send the result back */
#ifdef TCP_CORK
        rc = 1;
        setsockopt (data->sd, SOL_TCP, TCP_CORK, &rc, sizeof(rc));
#endif
        len = 0;
        while (len < sizeof(sf_cram_data_trvals)) { /* Header */
            rc = send (data->sd, (const void*)(((unsigned char*)data->trvals) + len),
                       sizeof(sf_cram_data_trvals) - len, MSG_NOSIGNAL);
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
        elen = data->trvals->n*data->nt*sizeof(float);
        if (data->kmah)
            elen *= 2;
        len = 0;
        while (len < elen) { /* Traces */
            rc = send (data->sd, (const void*)(((unsigned char*)traces) + len),
                       elen - len, MSG_NOSIGNAL);
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
    size_t i0, i1, n, nb;
    int ith = 1, nt, icpu, ncpu, tout;
    int port, nthreads, tmpfile = 0;
#ifdef LINUX
    int inet = 0;
#endif
    bool kmah;
    float trd, *traces = NULL;
    sf_file in, data = NULL, out;
    FILE *logf;
    char sbuffer[MAX_BUF];
    char *str = NULL;
    /* Fork variables */
    pid_t pid, sid;
    /* Server network variables */
    char *ip = NULL;
    int rc, on = 1, ijob, bsiz;
    int listen_sd = 0, new_sd, njobs = 0;
    struct sockaddr_in serv_addr, client_addr;
    struct timeval timeout;
    fd_set sset;
    socklen_t clen;
    /* Threading */
    sf_cramdd_work *rwork, *qjobs, *old_qjobs;
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

    /* Number of angular blocks to keep in memory per daemon */
    if (!sf_getint ("port", &port)) port = 18003;
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

    if (!sf_getint ("nthreads", &nthreads)) nthreads = 2*ncpu;
    /* Number of threads (connections) per daemon */
    if (!sf_getint ("timeout", &tout)) tout = 10;
    /* Inactivity time before shutdown (mins) */

    if (ith && 0 == (icpu % ith))
        sf_warning ("Running no more than %d threads on %s [CPU %d]", nthreads, ip, icpu);

    if (!sf_getstring ("data")) sf_error ("Need data=");
    /* Grid of supercells of local escape solutions */
    data = sf_input ("data");

    if (!sf_histbool (data, "KMAH", &kmah)) sf_error ("No KMAH= in data");
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in data");
    n = sf_leftsize (data, 1);
    if (kmah)
        n /= (size_t)2;
    /* Traces per daemon */
    trd = ith ? (float)n/((float)ncpu/(float)ith) : n;
    /* First and last trace indices to store */
    if (ith) {
        i0 = (int)(icpu/ith*trd);
        i1 = (int)((icpu/ith + 1)*trd);
        /* Add padding to compensate for possible truncation errors */
        nb = 1;
        if ((i1 - i0) > 100000)
            nb = 1000;
        else if ((i1 - i0) > 10000)
            nb = 100;
        else if ((i1 - i0) > 1000)
            nb = 10;
        i1 += nb;
        if (i1 >= n)
            i1 = n - 1;
        if (i0 != 0)
            i0 -= nb;
    } else {
        i0 = 0;
        i1 = n - 1;
    }

    sf_settype (out, SF_UCHAR);
    sf_putlargeint (out, "n1", sizeof (serv_addr) + 2*sizeof (size_t));
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Prestack data distributed access info");
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
    sf_putfloat (out, "trd", trd);
    sf_putstring (out, "Remote", ith ? "y" : "n");

    if (ith && 0 == (icpu % ith)) {
        sf_warning ("Trace range: %lu - %lu [CPU %d]", i0, i1, icpu);
    } else {
        i0 = n;
        i1 = 0;
    }

    sf_ucharwrite ((unsigned char*)&serv_addr, sizeof (serv_addr), out);
    sf_ucharwrite ((unsigned char*)&i0, sizeof (size_t), out);
    sf_ucharwrite ((unsigned char*)&i1, sizeof (size_t), out);

    if (ith && 0 == (icpu % ith)) {
        /* Temporary file for synchronizing forked processes */
        str = strdup ("/tmp/sfcramdd.XXXXXX");
        tmpfile = mkstemp (str);
        /* Daemonize */
        pid = fork ();
        if (pid < 0)
            sf_error ("fork() failed");
        if (0 == pid) {
            if (lockf (tmpfile, F_LOCK, 0) == -1)
                abort();
        }
        else
            sleep (1);
    } else
        pid = 1;

    /* Buffer seismic traces */
    if (ith && 0 == (icpu % ith) && 0 == pid) {
        /* Stream socket for incoming connections on */
        listen_sd = socket (AF_INET, SOCK_STREAM, 0);
        if (listen_sd < 0)
            sf_error ("socket() failed [CPU %d]", icpu);
        new_sd = connect (listen_sd, (struct sockaddr *)&serv_addr, sizeof(serv_addr));
        if (0 == new_sd) {
            sf_warning ("Daemon is already running [CPU %d]", icpu);
            close (new_sd);
            close (listen_sd);
            if (lockf (tmpfile, F_ULOCK, 1) == -1)
                abort();
            close(tmpfile);
            free (str);
            return 0;
        }
        n = i1 - i0 + 1;
        traces = sf_floatalloc (kmah ? n*(size_t)nt*(size_t)2 : n*(size_t)nt);
        sf_floatread (traces, kmah ? n*(size_t)nt*(size_t)2 : n*(size_t)nt, data);
        sf_warning ("%g Mb of traces buffered [CPU %d]",
                    kmah ? 1e-6*(float)n*(float)nt*2.0*sizeof(float)
                         : 1e-6*(float)n*(float)nt*sizeof(float), icpu);
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

    sf_fileclose (data);
    sf_fileclose (in);
    sf_fileclose (out);

    if (pid > 0)
        return 0;

    /* Change the file mode mask */
    umask (0);
    /* Create a new SID for the child process */
    sid = setsid ();
    if (sid < 0)
        sf_error ("setsid() failed [CPU %d]", icpu);
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
        sf_error ("setsockopt() failed [CPU %d]", icpu);
    }

    /* Set socket to be non-blocking; all of the sockets for
       the incoming connections will also be non-blocking */
    if (ioctl (listen_sd, FIONBIO, (char *)&on) < 0) {
        close (listen_sd);
        sf_error ("ioctl() failed [CPU %d]", icpu);
    }

    /* Bind the socket */
    memset (&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port = htons(port);
    if (bind (listen_sd,
              (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
        close (listen_sd);
        sf_error ("bind() failed [CPU %d]", icpu);
    }

    /* Set the listen back log */
    if (listen (listen_sd, 128) < 0) {
        close (listen_sd);
        sf_error ("listen() failed [CPU %d]", icpu);
    }

    pthread_mutex_init (&smutex, NULL);

    /* Open log file */
    snprintf (sbuffer, MAX_BUF, "sfcramdd_%s_%d.log", ip, icpu/ith);
    logf = fopen (sbuffer, "w+");
    fprintf (logf, "Listening on %s [CPU %d]\n", ip, icpu);
    fflush (logf);
    sf_warning ("Log file is %s [CPU %d]", sbuffer, icpu);

    /* Release the child before starting the server loop */
    if (lockf (tmpfile, F_ULOCK, 1) == -1)
        abort();
    close(tmpfile);
    free (str);
    fclose (stderr);

    /* Buffer for job requests */
    old_qjobs = NULL;
    qjobs = sf_cramdd_alloc_work (nthreads);
    /* Position in the buffer for the next client */
    ijob = 0;
    clen = sizeof(client_addr);

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
        /* Send buffer size */
        bsiz = kmah ? sizeof(sf_cram_data_trvals) + CRAM_TRBUF*nt*sizeof(float)*2
                    : sizeof(sf_cram_data_trvals) + CRAM_TRBUF*nt*sizeof(float);
        if (setsockopt (new_sd, SOL_SOCKET, SO_SNDBUF, &bsiz, sizeof(int)) < 0) {
            fprintf (logf, "Can not set SO_SNDBUF for a new connection\n");
            close (new_sd);
            break;
        }

        /* Connection with the new client has been established,
           prepare data for processing its requests */
        rwork = &qjobs[ijob];                    
        rwork->sd = new_sd;
        rwork->i0 = i0;
        rwork->i1 = i1;
        rwork->nt = nt;
        rwork->kmah = kmah;
        rwork->traces = traces;
        rwork->logf = logf;
        rwork->smutex = &smutex;
        rwork->njobs = &njobs;

        pthread_mutex_lock (&smutex);
        if (njobs < nthreads) {
            if (pthread_create (&pthread, NULL,
                                sf_cramdd_process_requests, (void*)rwork) != 0) {
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
                    sf_cramdd_free_work (old_qjobs);
                old_qjobs = qjobs;
                qjobs = sf_cramdd_alloc_work (nthreads);
                ijob = 0;
            }
            fprintf (logf, "Accepted client from %s, socket %d\n", ip, new_sd);
            fflush (logf);
        } else {
            fprintf (logf, "Reached maximum allowed number of threads, rejecting connection from %s\n", ip);
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
        sf_cramdd_free_work (old_qjobs);
    sf_cramdd_free_work (qjobs);
    free (traces);
    fclose (logf);

    return 0;
}

#else /* PTHREADS */

#include <rsf.h>

int main (int argc, char* argv[]) {
    fprintf (stderr, "sfcramdd can not work without pthreads; reconfigure Madagascar with pthreads support and recompile\n");
    return -1;
}
#endif

