/* Verify an SU dataset. */
/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

#include <stdio.h>

#include <rsf.h>

#include "verifytrace.h"

int main(int argc, char ** argv) 
{
    char * trialfile=NULL;
    char * compfile =NULL;
    int dim;
    float tol, rc, cref, amp, rdist;
    float e=0.0;
    float emax=0.0;
    int pf;
    segy trial;
    segy comp;
    FILE * fptrial=NULL;
    FILE * fpcomp=NULL;
    int ntr;

    sf_init(argc,argv);

    if (NULL == (trialfile = sf_getstring("trial")) || /*( trial trial file )*/
	NULL == (compfile  = sf_getstring("comp")))    /*( comp comparison file )*/
	sf_error("missing necessary inputs");

    if (!sf_getint("dim",&dim)) dim=3; /* dimension - leval values are 2 and 3 */
    if (!sf_getfloat("tol",&tol)) tol=0.05; /* tolerance for normalized max difference */
    if (!sf_getfloat("rc",&rc)) rc=0.1; /* nominal reflection coefficient, scales direct wave amp */
    if (!sf_getfloat("cref",&cref)) cref=1.5; /* nominal wave velocity */
    if (!sf_getfloat("amp",&amp)) amp=1.0; /* nominal amplitude of direct wave at reference distance */
    if (!sf_getfloat("rdist",&rdist)) rdist=1000.0; /* nominal reference distance */

    if (!(fptrial=fopen(trialfile,"r+")) ||
	!(fpcomp =fopen(compfile, "r+"))) 
	sf_error("failed to open files");
  
    ntr=0;
    emax=0.0;
    while (fgettr(fptrial,&trial) && fgettr(fpcomp,&comp)) {
	pf=verifytrace(&trial,&comp,dim,tol,rc,cref,amp,rdist,&e);
	sf_warning("trace %d: p/f = %d, scaled error = %e",ntr,pf,e);
	emax=fmax(e,emax);
	ntr++;
    }
    sf_warning("max scaled error over all traces = %e",emax);

    fclose(fptrial);
    fclose(fpcomp);

    exit(0);
}
