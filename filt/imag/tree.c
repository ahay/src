#include <math.h>

#include <rsf.h>

#include "tree.h"
#include "node.h"
#include "eno2.h"
#include "cell.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b))?(a):(b)
#endif

static Node Tree;
static NodeList Orphans;

static eno2 cvel;

static int nz, nx, na, nt, nax, naxz;
static float dz, dx, da, z0, x0, a0, **val;
static const float eps = 1.e-5;
static bool *accepted;

static void psnap (float* p, float* q, int* iq);

void tree_init (int order,
		int nz1, int nx1, int na1, int nt1, 
		float dz1, float dx1, float da1, 
		float z01, float x01, float a01, 
		float** vel, float** value) {
    nx = nx1; nz = nz1; na = na1; nt = nt1;
    dx = dx1; dz = dz1; da = da1;
    x0 = x01; z0 = z01; a0 = a01;

    nax = na*nx; naxz = nax*nz;

    cvel = eno2_init (order, nz, nx);
    eno2_set (cvel, vel); /* vel is slowness */

    val = value;
    accepted = sf_boolalloc(naxz);

    Orphans = CreateNodeList(nax);
    Tree = CreateNodes(naxz);
}

void tree_build(int method)
{
    int i, k, iz, ix, ia, kx, kz, ka, jx, jz;
    float x, z, p[2], a, v, v0, g0[2], g[2], s, sx, sz, t, *vk;
    bool onx, onz;
    Node node;

    for (kz=0; kz < nz; kz++) {
	sf_warning("Building %d of %d",kz+1,nz);
	for (kx=0; kx < nx; kx++) {
	
	    eno2_apply(cvel,kz,kx,0.,0.,&v0,g0,BOTH);
	    g0[1] /= dx;
	    g0[0] /= dz;

	    for (ka=0; ka < na; ka++) {
		k = ka + kx*na + kz*nax;
		node = Tree+k;

		a = a0 + ka*da;
		p[0] = -cos(a);
		p[1] = sin(a);

		/* get boundary conditions */
		if ((kx==0    && p[1] < 0.) ||
		    (kx==nx-1 && p[1] > 0.) ||
		    (kz==0    && p[0] < 0.) ||
		    (kz==nz-1 && p[0] > 0.)) {
		    AddNode(Orphans,node);
		    vk = val[k];
		    vk[0] = x0 + kx*dx;
		    vk[1] = z0 + kz*dz;
		    vk[2] = 0.;
		    vk[3] = cell_p2a (p);
		    accepted[k] = true;
		    continue;
		} else {
		    accepted[k] = false;
		}

		ia = ka;  
		x = 0.; ix=kx;
		z = 0.; iz=kz;		
		v = v0;
		g[0] = g0[0];
		g[1] = g0[1];

		switch (method) {
		    case 1:
			t = cell1_update2 (2, 0., v, p, g);
			/* p is normal vector now ||p|| = 1 */
	    
			cell1_intersect (g[1],x,dx/v,p[1],&sx,&jx);
			cell1_intersect (g[0],z,dz/v,p[0],&sz,&jz);
	    
			s = MIN(sx,sz);
	    
			t += cell1_update1 (2, s, v, p, g);
			/* p is slowness vector now ||p||=v */
	    
			if (s == sz) {
			    z = 0.; iz += jz;
			    x += p[1]*s/dx;
			} else {
			    x = 0.; ix += jx;
			    z += p[0]*s/dz;
			}
	    
			onz = cell_snap (&z,&iz,eps);
			onx = cell_snap (&x,&ix,eps);
	    
			eno2_apply(cvel,iz,ix,z,x,&v,g,BOTH);
			g[1] /= dx;
			g[0] /= dz;
	    
			t += cell1_update2 (2, s, v, p, g);
			/* p is normal vector now ||p||=1 */
			psnap (p,&a,&ia);
			break;
		    default:
			sf_error("Unknown method");
			break;
		}
  
		i = ia + ix*na + iz*nax;

		node->t = t;
		if (onz) { /* hits a z wall */
		    node->w1 = a;
		    node->w2 = x;
		    if (x != 1. && a != 1.) AddChild(Tree,i,0,node);
		    if (x != 1. && a != 0. && ia < na-1) 
			AddChild(Tree,i+1,1,node);
		    if (x != 0. && a != 1. && ix < nx-1) 
			AddChild(Tree,i+na,2,node);
		    if (x != 0. && a != 0. && ia < na-1 && ix < nx-1) 
			AddChild(Tree,i+na+1,3,node);
		} else { /* hits an x wall */
		    node->w1 = a;
		    node->w2 = z;
		    if (z != 1. && a != 1.) AddChild(Tree,i,0,node);
		    if (z != 1. && a != 0. && ia < na-1) 
			AddChild(Tree,i+1,1,node);
		    if (z != 0. && a != 1. && iz < nz-1) 
			AddChild(Tree,i+nax,2,node);
		    if (z != 0. && a != 0. && ia < na-1 && iz < nz-1) 
			AddChild(Tree,i+nax+1,3,node);
		}
	    }
	}
    }
}

void tree_print (void) {
    int k, ic;
    Node node;

    for (k=0; k < naxz; k++) {
	node = Tree+k;
	fprintf(stderr,"Node %d, nparents=%d, children: ",k,node->nparents);
	for (ic=0; ic < node->children->nitems; ic++) {
	    fprintf(stderr,"%d ",node->children->list[ic]-Tree);
	}
	fprintf(stderr,"\n");
    }
    fprintf(stderr,"Orphans: ");
    for (ic=0; ic < Orphans->nitems; ic++) {
	fprintf(stderr,"%d ",Orphans->list[ic]-Tree);
    }
    fprintf(stderr,"\n");
} 


void tree_traverse (void) {
    Node node, child;
    int k, i, j, n, nc, *parents;
    float w1, w2, *fk;

    for (n=0; n < Orphans->nitems; n++) {
	if (0==n%nax) fprintf(stderr,"Got %d of %d\n",n+1,naxz);

	node = Orphans->list[n];
	k = node - Tree;
	if (!accepted[k]) { /*evaluate node */      
	    fk = val[k];

	    fk[0] = 0.;
	    fk[1] = 0.;
	    fk[2] = node->t;
	    fk[3] = 0.;

	    w1 = node->w1;
	    w2 = node->w2;
	    parents = node->parents;

	    for (j=0; j < 4; j++) {
		if ((i=parents[0])>=0) fk[j] += (1.-w2)*(1.-w1)*val[i][j];
		if ((i=parents[1])>=0) fk[j] += (1.-w2)*    w1 *val[i][j];
		if ((i=parents[2])>=0) fk[j] +=     w2 *(1.-w1)*val[i][j];
		if ((i=parents[3])>=0) fk[j] +=     w2 *    w1 *val[i][j];
	    }
      
	    accepted[k] = true;
	}

	for (nc=0; nc < node->children->nitems; nc++) {
	    child = node->children->list[nc]; 
	    child->nparents--;
	    if (0==child->nparents) AddNode(Orphans,child);
	}
    }
}

void tree_close (void)
{
    FreeNodes(Tree,naxz);
    free(Orphans->list);
    free(Orphans);
}

static void psnap (float* p, float* q, int* iq) {
    int ia;
    float a2, a;

    a = cell_p2a(p);
    a2 = (a-a0)/da;
    ia = floor (a2); a2 -= ia;
    cell_snap (&a2, &ia, eps);

    if (ia < 0) {
	ia=0.; a2=0.; 
    } else if (ia > na-1 || (ia==na-1 && a2> 0.)) {
	ia=na-1; a2=0.;
    }

    a = a0+(ia+a2)*da;

    p[1] = sin(a);
    p[0] = -cos(a);

    *q = a2;
    *iq = ia;
}

