/* Tree structure for multiple arrivals. */
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

#include <math.h>
#include <assert.h>
#include <float.h>

#include <rsf.h>

#include "tree.h"
#include "node.h"

static Node Tree;
static NodeQueue Orphans;

static sf_eno2 cvel;

static int nz, nx, na, nax, naxz, nacc; /*, ii, jj; */
static float dz, dx, da, z0, x0, a0, **val;
static const float eps = 1.e-5;
static bool *accepted;

static void process_node (Node nd);
static void process_child (Node child);

void tree_init (int order     /* interpolation order */,
		int nz1       /* depth samples */, 
		int nx1       /* lateral samples */, 
		int na1       /* angle samples */, 
		float dz1     /* depth sampling */, 
		float dx1     /* lateral sampling */, 
		float da1     /* angle sampling */, 
		float z01     /* depth origin */, 
		float x01     /* lateral origin */,
		float a01     /* angle origin */,
		float** vel   /* slowness [nx][nz] */ , 
		float** value /* output [nx*nz*na][4] */) 
/*< initialize >*/
{
    int i;

    nx = nx1; nz = nz1; na = na1; 
    dx = dx1; dz = dz1; da = da1;
    x0 = x01; z0 = z01; a0 = a01;

    nax = na*nx; naxz = nax*nz;
    
    cvel = sf_eno2_init (order, nz, nx);
    sf_eno2_set (cvel, vel);

    val = value;
    accepted = sf_boolalloc(naxz);
    for (i=0; i < naxz; i++) {
	accepted[i] = false;
    }
    nacc = 0;

    Orphans = CreateNodeQueue();
    Tree = CreateNodes(naxz);
}

void tree_build(bool debug)
/*< Create a dependency tree >*/
{
    int i, k, iz, ix, ia, kx, kz, ka;
    float x, z, a, b, v, v0, p[2], g0[2], g[2], s, sx, sz, sp, t=0., *vk;
    bool onz, onp;
    Node node;

    for (kz=0; kz < nz; kz++) {
	sf_warning("Building %d of %d;",kz+1,nz);
	for (kx=0; kx < nx; kx++) {
	
	    sf_eno2_apply(cvel,kz,kx,0.,0.,&v0,g0,BOTH);
	    g0[1] /= dx;
	    g0[0] /= dz;

	    for (ka=0; ka < na; ka++) {
		k = ka + kx*na + kz*nax;
		node = Tree+k;

		a = a0+ka*da;
		p[0] = -cos(a);
		p[1] = sin(a);

		/* get boundary conditions */
		if ((kx==0    && p[1] <  FLT_EPSILON) ||
		    (kx==nx-1 && p[1] > -FLT_EPSILON) ||
		    (kz==0    && p[0] <  FLT_EPSILON) ||
		    (kz==nz-1 && p[0] > -FLT_EPSILON)) {
		    AddNode(Orphans,node);

		    vk = val[k];
		    vk[0] = x0 + kx*dx;
		    vk[1] = z0 + kz*dz;
		    vk[2] = 0.;
		    vk[3] = sf_cell_p2a(p);
		    accepted[k] = true;
		    nacc++;
		    continue;
		} else {
		    accepted[k] = false;
		}

		b = 0; ia = ka;  
		x = 0.; ix=kx;
		z = 0.; iz=kz;		
		v = v0;
		g[0] = g0[0];
		g[1] = g0[1];

		sx = v*p[1]/dx;
		sz = v*p[0]/dz;
		s = SF_MAX(fabsf(sx),fabsf(sz));

		sp = (g[0]*p[1]-g[1]*p[0])/da;
		s = SF_MAX(s,fabsf(sp));
	
		t = 0.5*v*v/s*(1.+(p[0]*g[0]+p[1]*g[1])/(3.*s));

		if (s == fabsf(sp)) {
		    b = 0.;
		    if (sp < 0) {
			ia--;
		    } else {
			ia++;
		    }
		    x += v*p[1]/(dx*s);
		    z += v*p[0]/(dz*s);
		} else if (s == fabsf(sz)) {
		    z = 0.; 
		    if (sz < 0.) {
			iz--;
		    } else {
			iz++;
		    }
		    x += v*p[1]/(dx*s);
		    b += (g[0]*p[1]-g[1]*p[0])/(da*s);
		} else {
		    x = 0.;
		    if (sx < 0.) {
			ix--;
		    } else {
			ix++;
		    }
		    z +=v*p[0]/(dz*s);
		    b += (g[0]*p[1]-g[1]*p[0])/(da*s);
		}
	    
		onz = sf_cell_snap (&z,&iz,eps);
		onp = sf_cell_snap (&b,&ia,eps);
		
		sf_eno2_apply(cvel,iz,ix,z,x,&v,g,BOTH);
		g[1] /= dx;
		g[0] /= dz;
	    
		a += (ia+b)*da;
		p[0] = -cos(a);
		p[1] = sin(a);

		t += 0.5*v*v/s*(1.-(p[0]*g[0]+p[1]*g[1])/(3.*s));
  		
		/* pathological exits */
		if (ix < 0 || ix > nx-1 ||
		    iz < 0 || iz > nz-1 ||
		    ia < 0 || ia > na-1) {
		    sf_warning("pathological exit (%d,%d,%d)",kz,kx,ka);

		    AddNode(Orphans,node);

		    vk = val[k];
		    vk[0] = x0 + kx*dx;
		    vk[1] = z0 + kz*dz;
		    vk[2] = 0.;
		    vk[3] = sf_cell_p2a(p);
		    accepted[k] = true;
		    nacc++;
		    continue;
		} 

		i = ia + ix*na + iz*nax;
		assert (i != k);
		
		node->t = t;
		if (onp) { /* hits a p wall */
		    node->w1 = z;
		    node->w2 = x;
		    if (x != 1. && z != 1.) 
			AddChild(Tree,i,0,0,node);
		    
		    if (iz == nz-1 || k == i+nax) {
			node->n1 = 1;
		    } else {
			node->n1 = 2;
			if (x != 1. && z != 0.)
			    AddChild(Tree,i+nax,0,1,node);		
		    }

		    if (ix == nx-1 || k == i+na) {
			node->n2 = 1;
		    } else {
			node->n2 = 2;
			if (x != 0. && z != 1.)
			    AddChild(Tree,i+na,1,0,node);
		    }

		    if (node->n1 == 2 && 
			node->n2 == 2 && 
			x != 0. && z != 0.) 
			AddChild(Tree,i+nax+na,1,1,node);
		} else if (onz) { /* hits a z wall */
		    node->w1 = b;
		    node->w2 = x;
		    if (x != 1. && b != 1.) 
			AddChild(Tree,i,0,0,node);
		    
		    if (ia == na-1 || k == i+1) {
			node->n1 = 1;
		    } else {
			node->n1 = 2;
			if (x != 1. && b != 0.)
			    AddChild(Tree,i+1,0,1,node);		
		    }

		    if (ix == nx-1 || k == i+na) {
			node->n2 = 1;
		    } else {
			node->n2 = 2;
			if (x != 0. && b != 1.)
			    AddChild(Tree,i+na,1,0,node);
		    }

		    if (node->n1 == 2 && 
			node->n2 == 2 && 
			x != 0. && b != 0.) 
			AddChild(Tree,i+na+1,1,1,node);
		} else { /* hits an x wall */
		    node->w1 = b;
		    node->w2 = z;
		    if (z != 1. && b != 1.)
			AddChild(Tree,i,0,0,node);

		    if (ia == na-1 || k == i+1) {
			node->n1 = 1;
		    } else {
			node->n1 = 2;
			if (z != 1. && b != 0.) 
			    AddChild(Tree,i+1,0,1,node);
		    }

		    if (iz == nz-1 || k == i+nax) {
			node->n2 = 1;
		    } else {
			node->n2 = 2;
			if (z != 0. && b != 1.) 
			    AddChild(Tree,i+nax,1,0,node);
		    }

		    if (node->n1 == 2 && node->n2 == 2 && 
			z != 0. && b != 0.) 
			AddChild(Tree,i+nax+1,1,1,node);
		}
	    }
	}
    } 
    sf_warning(".");

    if (debug) tree_print();

    TraverseDeleteQueue (Orphans,process_node);

    if (debug) tree_print();
	
    if (nacc == naxz) return;

    sf_warning("Found %d < %d, entering cycle resolution",nacc,naxz);
	
    /* 
    sf_pqueue_init (naxz);
    sf_pqueue_start ();

    check_front();

    while (nacc < naxz) {
	vk = sf_pqueue_extract ();
	if (NULL == vk) {
	    sf_warning("Heap exhausted: %d accepted",nacc);
	    tree_print();
	    return;
	}
	node = Tree + ((vk-val[0]-2)/4);
	TraverseDeleteQueue(node->children,orphanize);
    }
    */
}

/*
static void catch_node (Node node) {
    if (ii == node-Tree) sf_error("got it %d!",jj);
    jj++;
    TraverseQueue(node->children,catch_node);
}
*/

static void print_node (Node node) {
    int k, kx, kz, ka;

    k = node-Tree;
    ka = k;
    kz = ka/nax; ka -= nax*kz;
    kx = ka/na;  ka -= na*kx;

    if (accepted[k]) {
	fprintf(stderr,"[%d %d %d] ",ka+1,kx+1,kz+1);
    } else {
	fprintf(stderr,"{%d %d %d} ",ka+1,kx+1,kz+1);
    }
}

void tree_print (void) 
/*< Print out the tree (for debugging >*/
{
    int k;
    Node node;

    for (k=0; k < naxz; k++) {
	node = Tree+k;
	fprintf(stderr,"Node ");
	print_node(node);
	fprintf(stderr,"nparents=%d, children: ",node->nparents);
	if (NULL != node->children)
	    TraverseQueue (node->children,print_node);
	fprintf(stderr,"\n");
    }
    fprintf(stderr,"Orphans: ");
    TraverseQueue (Orphans,print_node);
    fprintf(stderr,"\n");
} 

static void process_node (Node nd) {
    static int n=0;
    int k, i, j, k1, k2, **parents;
    float x, w1[2], w2[2], *fk;

    if (0==n%nax) sf_warning("Got %d of %d",n+1,naxz);
    n++;

    k = nd - Tree;
    if (!accepted[k]) { /*evaluate node */      
	fk = val[k];

	fk[0] = 0.;
	fk[1] = 0.;
	fk[2] = nd->t;
	fk[3] = 0.;
	
	parents = nd->parents;
	
	x = nd->w1; 
	switch (nd->n1) {
	    case 1:
		w1[0] = 1;
		break;
	    case 2:
		w1[0] = 1.-x; 
		w1[1] = x;
		break;
	}
	
	x = nd->w2; 
	switch (nd->n2) {
	    case 1:
		w2[0] = 1;
		break;
	    case 2:
		w2[0] = 1.-x; 
		w2[1] = x;
		break;
	}
	
	for (k2=0; k2 < nd->n2; k2++) {		  
	    for (k1=0; k1 < nd->n1; k1++) {
		i = parents[k2][k1];
		if (i >= 0) {					
		    x = w2[k2]*w1[k1];
		    
		    for (j=0; j < 4; j++) {
			fk[j] += x*val[i][j];
		    }
		}
	    }
	}
	
	accepted[k] = true;
	nacc++;
    }
    
    free (nd->parents);
    TraverseQueue (nd->children,process_child);
}

static void process_child (Node child) {
    child->nparents--;
    if (0==child->nparents) AddNode(Orphans,child);
}

void tree_close (void)
/*< Free allocated storage >*/
{
    FreeNodes(Tree,naxz);
}

/* 	$Id: tree.c 1575 2005-11-21 14:09:06Z fomels $	 */
