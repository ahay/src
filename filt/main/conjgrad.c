/* Generic conjugate-gradient solver for linear inversion */
/*
  Copyright (C) 2005 University of Texas at Austin
  
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
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
/*^*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, iter, niter, p[6][2];
    char *name;
    float *buf, *buf2;
    double rn, rnp, alpha, beta;
    pid_t pid[5]={1,1,1,1,1};
    off_t nm, nd, msiz, dsiz, pos;
    size_t nbuf, mbuf, dbuf;
    FILE *x, *rr, *g, *s;
    sf_file mod, dat, out, from, to;

    sf_init(argc,argv);
    dat = sf_input("in");
    mod = sf_input("mod");

    if (SF_FLOAT != sf_gettype(mod) ||
	SF_FLOAT != sf_gettype(dat)) 
	sf_error("Need float type in mod and dat");
  
    free(argv[0]);
    for (i=0; i < argc-1; i++) {
	argv[i]=argv[i+1];
    }
    argv[argc-1] = sf_charalloc(6);
    snprintf(argv[argc-1],6,"adj=X");

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of iterations */

    rr = sf_tempfile(&name,"w+b");
    x  = sf_tempfile(&name,"w+b");
    g  = sf_tempfile(&name,"w+b");
    s  = sf_tempfile(&name,"w+b");

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    /* I/O buffers */
    nbuf = BUFSIZ/sizeof(float);
    buf  = sf_floatalloc(nbuf);
    buf2 = sf_floatalloc(nbuf);

    /* r=-dat */

    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	if (dsiz < dbuf) dbuf=dsiz;

	sf_floatread(buf,dbuf,dat);
	for (i=0; i < dbuf; i++) {
	    buf[i]=-buf[i];
	}

	fwrite(buf,sizeof(float),dbuf,rr);
    }
    
    /* x=0 */

    for (i=0; i < nbuf; i++) {
	buf[i]=0.;
    }
    
    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	if (msiz < mbuf) mbuf=msiz;
	fwrite(buf,sizeof(float),mbuf,x);
    }

    for (i=0; i < 6; i++) { /* make six pipes */
	if (pipe(p[i]) < 0) sf_error("pipe error:");
    }

    for (iter=0; iter < niter; iter++) {
	sf_warning("iter %d of %d",iter+1,niter);
	

	for (i=0; i < 5; i++) { /* fork five children */
	    if ((pid[i] = fork()) < 0) sf_error("fork error:");
	    if (0 == pid[i]) break;
	}
	
	if (0 == pid[0]) {	
	    /* feeds rr to p[0] */

	    sf_warning("child 0 start");

	    close(p[0][0]);
	    close(STDOUT_FILENO);
	    dup(p[0][1]);

	    to = sf_output("out");

	    if (0 > fseeko(rr,0,SEEK_SET))
		sf_error ("seek problem");

	    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
		if (dsiz < dbuf) dbuf=dsiz;

		if (dbuf != fread(buf,sizeof(float),dbuf,rr)) sf_error("read error");

		sf_floatwrite(buf,dbuf,to);
	    }

	    sf_warning("child 0 end");

	    _exit(0);
	}
 
	
	if (0==pid[1]) {
	    /* reads from p[0], runs the adjoint, and writes to p[1] */

	    sf_warning("child 1 start");

	    close(p[0][1]);
	    close(STDIN_FILENO);
	    dup(p[0][0]);
	    
	    close(p[1][0]);
	    close(STDOUT_FILENO);
	    dup(p[1][1]);

	    sf_warning("child 1 launch");
	    
	    argv[argc-1][4]='1';
	    execvp(argv[0],argv);
	    _exit(1);
	}


	if (0==pid[2]) {
	    /* reads gradient from p[1], uses it, and writes to p[2] */
	
	    sf_warning("child 2 start");

	    close(p[1][1]);
	    close(STDIN_FILENO);
	    dup(p[1][0]);
	    from = sf_input("in");
	    
	    close(p[2][0]);
	    close(STDOUT_FILENO);
	    dup(p[2][1]);
	    to = sf_output("out");
	    sf_fileflush(to,mod);

	    if (0 > fseeko(g,0,SEEK_SET))
		sf_error ("seek problem");

	    /* rn = ddot(g) */
	    
	    rn = 0.;

	    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
		if (msiz < mbuf) mbuf=msiz;
		sf_floatread(buf,mbuf,from);

		for (i=0; i < mbuf; i++) {
		    rn += (double) buf[i] * buf[i];
		}

		fwrite(buf,sizeof(float),mbuf,g);
		sf_floatwrite(buf,mbuf,to);
	    }

	    if (0==iter) {
		alpha = 0.;
	    } else {
		if (0 > fseeko(s,0,SEEK_SET))
		    sf_error ("seek problem");

		if (1 != fread(&rnp,sizeof(float),1,s)) sf_error("read error");

		alpha = rn/rnp;

		if (0 > fseeko(s,0,SEEK_SET))
		    sf_error ("seek problem");
	    }

	    fwrite(&rn,sizeof(float),1,s);
	    fflush(s);

	    write(p[4][1],&alpha,sizeof(float));
	    
	    if (0 > fseeko(g,0,SEEK_SET))
		sf_error ("seek problem");

	    /* s = g + alpha*s */
	    
	    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
		if (msiz < mbuf) mbuf=msiz;

		if (mbuf != fread(buf,sizeof(float),mbuf,g)) sf_error("read error");

		if (iter > 0) {
		    pos = ftello(s);

		    if (mbuf != fread(buf2,sizeof(float),mbuf,s))  sf_error("read error");

		    if (0 > fseeko(s,pos,SEEK_SET))
			sf_error ("seek problem");

		    for (i=0; i < mbuf; i++) {
			buf[i] += alpha * buf2[i];
		    }
		} 

		fwrite(buf,sizeof(float),mbuf,s);
	    }

	    sf_warning("child 2 done");

	    _exit(2);
	}

	if (0==pid[3]) {
	    /* reads from p[2], runs the forward, and writes to p[3] */

	    sf_warning("child 3 start");

	    close(p[2][1]);
	    close(STDIN_FILENO);
	    dup(p[2][0]);
	    
	    close(p[3][0]);
	    close(STDOUT_FILENO);
	    dup(p[3][1]);

	    sf_warning("child 3 launch");
	    
	    argv[argc-1][4]='0';
	    execvp(argv[0],argv);
	    _exit(3);
	}

	if (0==pid[4]) {
	    /* reads conj. grad from p[3], does computations with it */

	    sf_warning("child 4 start");

	    close(p[3][1]);
	    close(STDIN_FILENO);
	    dup(p[3][0]);
	    from = sf_input("in");
	
	    read(p[4][0],&alpha,sizeof(float));

	    if (0 > fseeko(s,0,SEEK_SET))
		sf_error ("seek problem");

	    if (1 != fread(&rn,sizeof(float),1,s)) sf_error("read error");

	    if (0 > fseeko(s,sizeof(float)*nm,SEEK_CUR))
		sf_error ("seek problem");

	    /* beta = ddot(ss) */

	    beta = 0.;

	    /* ss = gg + alpha * ss */

	    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
		if (dsiz < dbuf) dbuf=dsiz;
		sf_floatread(buf,dbuf,from);

		if (iter > 0) {
		    pos = ftello(s);

		    if (dbuf != fread(buf2,sizeof(float),dbuf,s)) sf_error("read error");

		    if (0 > fseeko(s,pos,SEEK_SET))
			sf_error ("seek problem");
		    
		    for (i=0; i < dbuf; i++) {
			buf[i] += alpha * buf2[i];
		    }
		}

		for (i=0; i < dbuf; i++) {
		    beta += (double) buf[i] * buf[i];
		}
		
		fwrite(buf,sizeof(float),dbuf,s);
	    }

	    alpha = - rn/beta;

	    write(p[5][1],&alpha,sizeof(float));

	    sf_warning("child 4 done");

	    _exit(4);
	}

	for (i=0; i < 5; i++) { 
	    if (0 == pid[i]) sf_error("A child is alive");
	}    

	/* parent updates x and r */
	read(p[5][0],&alpha,sizeof(float));

	if (0 > fseeko(s,sizeof(float),SEEK_SET))
	    sf_error ("seek problem");

	if (0 > fseeko(x,0,SEEK_SET))
		sf_error ("seek problem");

	for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	    if (msiz < mbuf) mbuf=msiz;

	    pos = ftello(x);

	    if (mbuf != fread(buf,sizeof(float),mbuf,x) ||
		mbuf != fread(buf2,sizeof(float),mbuf,s)) sf_error("read error (x)");

	    for (i=0; i < mbuf; i++) {
		buf[i] += alpha * buf2[i];
	    }

	    if (0 > fseeko(x,pos,SEEK_SET))
		sf_error ("seek problem");
	    
	    fwrite(buf,sizeof(float),mbuf,x);
	}

	if (0 > fseeko(rr,0,SEEK_SET))
		sf_error ("seek problem");

	for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	    if (dsiz < dbuf) dbuf=dsiz;

	    pos = ftello(rr);

	    if (dbuf != fread(buf,sizeof(float),dbuf,rr) ||
		dbuf != fread(buf2,sizeof(float),dbuf,s)) sf_error("read error (rr)");

	    for (i=0; i < dbuf; i++) {
		buf[i] += alpha * buf2[i];
	    }

	    if (0 > fseeko(rr,pos,SEEK_SET))
		sf_error ("seek problem");

	    
	    fwrite(buf,sizeof(float),dbuf,rr);
	}
    }

    /* write x to out */  
    out = sf_output("out");
    sf_fileflush(out,mod);
  
    if (0 > fseeko(x,0,SEEK_SET))
	sf_error ("seek problem");

    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	if (msiz < mbuf) mbuf=msiz;

	if (mbuf != fread(buf,sizeof(float),mbuf,x)) sf_error("read error (last)");
	
	sf_floatwrite(buf,mbuf,out);
    }

    exit(0);
}



