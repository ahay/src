/* Generic dot-product test for linear operators with adjoints */
/*
  Copyright (C) 2004 University of Texas at Austin

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
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>

#include <rsf.h>

#define DUP(a) if (dup(a) < 0) sf_error("dup error:")

int main(int argc, char* argv[])
{
    int p[4][2], i, im, id, status;
    unsigned long mseed, dseed;
    off_t nm, nd, msiz, dsiz;
    size_t nbuf, mbuf, dbuf;
    float *buf;
    double dp;
    pid_t pid[6]={1,1,1,1,1,1};
    sf_file mod=NULL;
    sf_file dat=NULL;
    sf_file pip=NULL;

    sf_init(argc,argv);

    mod = sf_input("mod");
    dat = sf_input("dat");

    if (SF_FLOAT != sf_gettype(mod) ||
	SF_FLOAT != sf_gettype(dat))
	sf_error("Need float type in mod and dat");

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    nbuf = BUFSIZ/sizeof(float);
    buf = sf_floatalloc(nbuf);

    mseed = (unsigned long) time(NULL);
    init_genrand(mseed);
    mseed = genrand_int32();
    dseed = genrand_int32();

    for (i=0; i < argc-1; i++) {
	argv[i]=argv[i+1];
    }
    argv[argc-1] = sf_charalloc(6);
    snprintf(argv[argc-1],6,"adj=X");

    for (i=0; i < 4; i++) { /* make four pipes */
	if (pipe(p[i]) < 0) sf_error("pipe error:");
    }

    for (i=0; i < 6; i++) { /* fork six children */
	if ((pid[i] = fork()) < 0) sf_error("fork error:");
	if (0 == pid[i]) break;
    }

    if (0 == pid[0]) {	
	/* makes random model and writes it to p[0] */

	close(p[0][0]);
	close(STDOUT_FILENO);
	DUP(p[0][1]);

	pip = sf_output("out");
	sf_fileflush(pip,mod);

	init_genrand(mseed);
	for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	    if (msiz < mbuf) mbuf=msiz;

	    sf_random(mbuf,buf);

	    sf_floatwrite(buf,mbuf,pip);
	}
    } 

    if (0 == pid[1]) {
	/* reads from p[0], runs the program, and writes to p[1] */

	close(p[0][1]);
	close(STDIN_FILENO);
	DUP(p[0][0]);

	close(p[1][0]);
	close(STDOUT_FILENO);
	DUP(p[1][1]);

	argv[argc-1][4]='0';
	execvp(argv[0],argv);

	_exit(1);
    }

    if (0 == pid[2]) {
	/* reads from p[1] and multiplies it with random data */
	
	close(p[1][1]);
	close(STDIN_FILENO);
	DUP(p[1][0]);

	pip = sf_input("in");

	init_genrand(dseed);
	dp = 0.;
	for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	    if (dsiz < dbuf) dbuf=dsiz;

	    sf_floatread(buf,dbuf,pip);
	    for (id=0; id < dbuf; id++) {
		dp += buf[id]*genrand_real1 ();
	    }	
	}
	sf_warning(" L[m]*d=%g",dp);

	_exit(2);
    }

    if (0 == pid[3]) {	
	/* makes random data and writes it to p[2] */

	close(p[2][0]);
	close(STDOUT_FILENO);
	DUP(p[2][1]);

	pip = sf_output("out");
	sf_fileflush(pip,dat);

	init_genrand(dseed);
	for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	    if (dsiz < dbuf) dbuf=dsiz;

	    sf_random(dbuf,buf);

	    sf_floatwrite(buf,dbuf,pip);
	}
    } 

    if (0 == pid[4]) {
	/* reads from p[2], runs the adjoint, and writes to p[3] */

	close(p[2][1]);
	close(STDIN_FILENO);
	DUP(p[2][0]);

	close(p[3][0]);
	close(STDOUT_FILENO);
	DUP(p[3][1]);

	argv[argc-1][4]='1';
	execvp(argv[0],argv);

	_exit(4);
    }

    if (0 == pid[5]) {
	/* reads from p[3] and multiplies it with random model */
	
	close(p[3][1]);
	close(STDIN_FILENO);
	DUP(p[3][0]);

	pip = sf_input("in");

	init_genrand(mseed);
	dp = 0.;
	for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	    if (msiz < mbuf) mbuf=msiz;

	    sf_floatread(buf,mbuf,pip);


	    for (im=0; im < mbuf; im++) {
		dp += buf[im]*genrand_real1 ();
	    }	
	}
	sf_warning("L'[d]*m=%g",dp);
	
	_exit(5);
    }

    for (i=0; i < 6; i++) {
	if (0 == pid[i]) break;
    }

    if (6==i) {
	/* parent waits */
	waitpid(pid[2],&status,0);
	waitpid(pid[5],&status,0);
    
	exit(0);
    }
}
