/* OpenMP wrapper for embarassingly parallel jobs. */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include <omp.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int axis, axis2, rank, nodes, ndim, jobs;
    off_t n[SF_MAX_DIM];
    char *iname=NULL, **cmdline;
    FILE *tmp;
    sf_file inp, out, inp2;

#pragma omp parallel
    {
	nodes = omp_get_num_threads();
	if (1 >= nodes) 
	    nodes = omp_get_num_procs(); 
    }

    /* master node */
    sf_init(argc,argv);

    inp = sf_input("in");
    out = sf_output("out");

    ndim = sf_largefiledims (inp,n);
    if (!sf_getint("split",&axis)) axis=ndim;
    /* axis to split */

    tmp = sf_tempfile(&iname,"w+b");
    fclose(tmp);

    inp2 = sf_output(iname);
    sf_cp(inp,inp2);
    sf_fileclose(inp2);

    inp2 = sf_input(iname);

    cmdline = sf_split(inp2,axis,nodes+1,&jobs,ndim,n,argc,argv);  
    sf_warning("Running %d threads",jobs);

#pragma omp parallel
    {
	omp_set_num_threads(jobs);
    }

#pragma omp parallel private(rank) shared(cmdline)
    {
	rank = omp_get_thread_num();
	if (rank < jobs) {
	    fprintf(stderr,"CPU %d: %s\n",rank,cmdline[rank]); 
	    sf_system(cmdline[rank]);
	}
    }
    
    if (!sf_getint("join",&axis2)) axis2=axis;
    /* axis to join (0 means add) */
    
    sf_out(out,jobs,axis2,iname);
    sf_rm(iname,true,false,false);

    if (axis2 > 0) {
	sf_join(out,axis2,jobs);
    } else {
	sf_add(out,jobs);
    }

    exit(0);
}
