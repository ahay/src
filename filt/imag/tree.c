#include <math.h>
#include <assert.h>

#include <rsf.h>

#include "tree.h"
#include "node.h"
#include "eno2.h"
#include "cell.h"
#include "pqueue.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b))?(a):(b)
#endif

static Node Tree;
static NodeQueue Orphans;

static eno2 cvel;

static int nz, nx, na, nt, nax, naxz, order, nacc, ii, jj;
static float dz, dx, da, z0, x0, a0, **val;
static const float eps = 1.e-5;
static bool *accepted;

static void psnap (float* p, float* q, int* iq);

static void process_node (Node nd);
static void process_child (Node child);
static void check_front (void);

void tree_init (int order1,
		int nz1, int nx1, int na1, int nt1, 
		float dz1, float dx1, float da1, 
		float z01, float x01, float a01, 
		float** vel, float** value) {
    nx = nx1; nz = nz1; na = na1; nt = nt1;
    dx = dx1; dz = dz1; da = da1;
    x0 = x01; z0 = z01; a0 = a01;

    nax = na*nx; naxz = nax*nz;
    order = order1;
    
    cvel = eno2_init (order, nz, nx);
    eno2_set (cvel, vel); /* Del is slowness */

    val = value;
    accepted = sf_boolalloc(naxz);
    nacc = 0;

    Orphans = CreateNodeQueue();
    Tree = CreateNodes(naxz,order);
}

void tree_build(void)
{
    int i, k, iz, ix, ia, kx, kz, ka, jx, jz;
    float x, z, p[2], a, v, v0, g0[2], g[2], s, sx, sz, t=0., *vk;
    bool onx, onz=false;
    Node node;

    sf_warning("Method=%d",order);
    for (kz=0; kz < nz; kz++) {
	sf_warning("Building %d of %d",kz+1,nz);
	for (kx=0; kx < nx; kx++) {
	
	    eno2_apply(cvel,kz,kx,0.,0.,&v0,g0,BOTH);
	    g0[1] /= dx;
	    g0[0] /= dz;

	    for (ka=0; ka < na; ka++) {
		k = ka + kx*na + kz*nax;
		node = Tree+k;

		/*** debug ***
		sf_warning("node %d %d %d",ka,kx,kz);
		*************/

		a = a0 + ka*da;
		p[0] = -cos(a);
		p[1] = sin(a);

		/* get boundary conditions */
		if ((kx==0    && p[1] < 0.) ||
		    (kx==nx-1 && p[1] > 0.) ||
		    (kz==0    && p[0] < 0.) ||
		    (kz==nz-1 && p[0] > 0.)) {
		    AddNode(Orphans,node);

		    /*** debug ***
		    sf_warning("orphan");
		    *************/

		    vk = val[k];
		    vk[0] = x0 + kx*dx;
		    vk[1] = z0 + kz*dz;
		    vk[2] = 0.;
		    vk[3] = cell_p2a (p);
		    accepted[k] = true;
		    nacc++;
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

		switch (order) {
		    case 2:
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
		    case 3:
			t = cell_update2 (2, 0., v, p, g);
			/* p is normal vector now ||p|| = 1 */
	    
			cell_intersect (g[1],x,dx/v,p[1],&sx,&jx);
			cell_intersect (g[0],z,dz/v,p[0],&sz,&jz);
	    
			s = MIN(sx,sz);
	    
			t += cell_update1 (2, s, v, p, g);
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
	    
			t += cell_update2 (2, s, v, p, g);
			/* p is normal vector now ||p||=1 */
			psnap (p,&a,&ia);
			break;
		    default:
			sf_error("Unknown method");
			break;
		}
  		
		/* pathological exits */
		if (ix < 0 || ix > nx-1 ||
		    iz < 0 || iz > nz-1) {
		    AddNode(Orphans,node);

		    /*** debug ***
		    sf_warning("orphan 2");
		    *************/

		    vk = val[k];
		    vk[0] = x0 + kx*dx;
		    vk[1] = z0 + kz*dz;
		    vk[2] = 0.;
		    vk[3] = cell_p2a (p);
		    accepted[k] = true;
		    nacc++;
		    continue;
		} 

		i = ia + ix*na + iz*nax;
		assert (i != k);

		node->t = t;
		if (onz) { /* hits a z wall */
		    node->w1 = a;
		    node->w2 = x;
		    if (x != 1. && a != 1.) 
			AddChild(Tree,i,0,0,node);
		    
		    if (ia == na-1 || k == i+1) {
			node->n1 = 1;
		    } else {
			node->n1 = order;
			if (x != 1. && a != 0.)
			    AddChild(Tree,i+1,0,1,node);		
		    }

		    if (ix == nx-1 || k == i+na) {
			node->n2 = 1;
		    } else {
			node->n2 = order;
			if (x != 0. && a != 1.)
			    AddChild(Tree,i+na,1,0,node);
		    }

		    if (node->n1 == order && 
			node->n2 == order && 
			x != 0. && a != 0.) 
			AddChild(Tree,i+na+1,1,1,node);

		    if (3==order) {
			if (ia == 0 || k == i-1) {
			    if (node->n1 > 1) node->n1 = 2;
			} else if (x != 1. && a != 0. && a != 1.) {
				AddChild(Tree,i-1,0,2,node);
			}

			if (ix == 0 || k == i-na) {
			    if (node->n2 > 1) node->n2 = 2;
			} else if (x != 0. && x != 1. && a != 1.) {
			    AddChild(Tree,i-na,2,0,node);
			}

			if (x != 0. && a != 0. && 
			    node->n1 > 1 && node->n2 > 1) {
			    if (node->n1 == 3 && a != 1.) {
				if (k == i+na-1) {
				    node->n1 = 2;
				} else {
				    AddChild(Tree,i+na-1,1,2,node);

				    if (node->n2 == 3 && x != 1.) {
					if (k == i-na-1) {
					    node->n2 = 2;
					} else {
					    AddChild(Tree,i-na-1,2,2,node);
					}
				    }
				}
			    }
			    if (node->n2 == 3 && x != 1.) {
				if (k == i-na+1) {
				    node->n2 = 2;
				} else {
				    AddChild(Tree,i-na+1,2,1,node);
				}
			    }
			}			
		    }
		} else { /* hits an x wall */
		    node->w1 = a;
		    node->w2 = z;
		    if (z != 1. && a != 1.)
			AddChild(Tree,i,0,0,node);

		    if (ia == na-1 || k == i+1) {
			node->n1 = 1;
		    } else {
			node->n1 = order;
			if (z != 1. && a != 0.) 
			    AddChild(Tree,i+1,0,1,node);
		    }

		    if (iz == nz-1 || k == i+nax) {
			node->n2 = 1;
		    } else {
			node->n2 = order;
			if (z != 0. && a != 1.) 
			    AddChild(Tree,i+nax,1,0,node);
		    }

		    if (node->n1 == order && node->n2 == order && 
			z != 0. && a != 0.) 
			AddChild(Tree,i+nax+1,1,1,node);


		    if (3==order) {
			if (ia == 0 || k == i-1) {
			    if (node->n1 > 1) node->n1 = 2;
			} else if (z != 1. && a != 0. && a != 1.) {
			    AddChild(Tree,i-1,0,2,node);
			}

			if (iz == 0 || k == i-nax) {
			    if (node->n2 > 1) node->n2 = 2;
			} else if (z != 0. && z != 1. && a != 1.) {
			    AddChild(Tree,i-nax,2,0,node);
			}
			
			if (z != 0. && a != 0. && 
			    node->n1 > 1 && node->n2 > 1) {
			    if (node->n1 == 3 && a != 1.) { 
				if (k == i+nax-1) {
				    node->n1 = 2;
				} else {
				    AddChild(Tree,i+nax-1,1,2,node);

				    if (node->n2 == 3 && z != 1.) {
					if (k == i-nax-1) {
					    node->n2 = 2;
					} else {
					    AddChild(Tree,i-nax-1,2,2,node);
					}
				    }
				}
			    }
			    if (node->n2 == 3 && z != 1.) { 
				if (k == i-nax+1) {
				    node->n2 = 2;
				} else {
				    AddChild(Tree,i-nax+1,2,1,node);
				}
			    }
			}
		    }
		}
	    }
	}
    }

/*    tree_print(); */

    TraverseQueue (Orphans,process_node);
	
    if (nacc == naxz) return;

    sf_warning("Found %d < %d, entering cycle resolution",nacc,naxz);

    FreeNodeQueue (Orphans);
    Orphans = CreateNodeQueue();
	
    pqueue_init (naxz-nacc);
    pqueue_start ();

    check_front();
}

static void check_front (void) {
    int k;
    bool atfront=false;
    NodeCell cell;
    Node node;

    for (k=0; k < naxz; k++) {
	if (accepted[k]) {
	    node = Tree+k;
	    for (cell = node->children->head; 
		 NULL != cell; 
		 cell = cell->link) {
		atfront = (cell->node->nparents > 0);
		if (atfront) break;
	    }
	    if (atfront) {
		pqueue_insert(val[k]+2);
		sf_warning("%d",k);
	    }
	}
    }
}

static void catch_node (Node node) {
    if (ii == node-Tree) sf_error("got it %d!",jj);
    jj++;
    TraverseQueue(node->children,catch_node);
}

static void catch(int k) {
    Node node;

    ii = k;
    jj = 0;
    node = Tree+k;
    TraverseQueue(node->children,catch_node);
}

static void print_node (Node node) {
    int k, kx, kz, ka;
    
    k = node-Tree;
    kz = k/nax; k -= nax*kz;
    kx = k/na;  k -= na*kx;
    ka = k;

    fprintf(stderr,"[%d %d %d] ",ka+1,kx+1,kz+1);
}

/* print_queue */
void tree_print (void) {
    int k;
    Node node;

    catch(0);

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
    float x, w1[3], w2[3], *fk;

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
	    case 3:
		w1[0] = 1.-x*x; 
		w1[1] = 0.5*x*(x+1.); 
		w1[2] = 0.5*x*(x-1.);
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
	    case 3:
		w2[0] = 1.-x*x; 
		w2[1] = 0.5*x*(x+1.); 
		w2[2] = 0.5*x*(x-1.);
		break;
	}
	
	/**** debug ****
	      kz = k/nax; 
	      kx = (k-kz*nax)/na;
	      ka = k-kz*nax-kx*na;
	      sf_warning("node %d %d %d",ka,kx,kz);
	***************/
	

	for (k2=0; k2 < nd->n2; k2++) {		  
	    for (k1=0; k1 < nd->n1; k1++) {
		i = parents[k2][k1];
		if (i >= 0) {					
		    x = w2[k2]*w1[k1];
		    
		    /**** debug ****
			  iz = i/nax; 
			  ix = (i-iz*nax)/na;
			  ia = i-iz*nax-ix*na;
			  sf_warning("weight %d %d %d: %g",ia,ix,iz,x);
		    ***************/
		    
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
    FreeNodeQueue (nd->children);
    nd->children = NULL;
}

static void process_child (Node child) {
    child->nparents--;
    if (0==child->nparents) AddNode(Orphans,child);
}

void tree_close (void)
{
    FreeNodes(Tree,naxz);
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

