/* Surface-consistent decomposition */
#include <rsf.h>

int main(int argc, char* argv[])
{
    bool adj, verb;
    int nd, nm, nx, im, min, max, id, i, ix, sx;
    int **indx, *size;
    float *model, *data;
    sf_file inp, index, out;

    sf_init(argc,argv);
  
    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    inp = sf_input("in");
    if (SF_FLOAT != sf_gettype(inp)) 
	sf_error("Need float input");

    out = sf_output("out");

    index = sf_input("index");
    if (SF_INT != sf_gettype(index)) 
	sf_error("Need int index");

    if (!sf_histint(index,"n1",&nd)) 
	sf_error("No n1= in index");
    nm = sf_leftsize(index,1);

    if (adj) {
	if (nd != sf_filesize(inp)) 
	    sf_error("Wrong data size");
    } else {
	sf_putint(out,"n1",nd);
    }

    data = sf_floatalloc(nd);
    indx = sf_intalloc2(nd,nm);
    size = sf_intalloc(nm);

    sf_intread(indx[0],nd*nm,index);

    nx = 0;
    for (im=0; im < nm; im++) {
	min = max = indx[im][0];
	for (id=1; id < nd; id++) {
	    i = indx[im][id];
	    if (i < min) min=i;
	    if (i > max) max=i;
	}
	if (min) {
	    for (id=0; id < nd; id++) {
		indx[im][id] -= min;
	    }
	}
	size[im]=max-min+1;
	nx += size[im];
	if (verb) sf_warning("size%d=%d",im+1,size[im]);
    }

    if (adj) {
	sf_putint(out,"n1",nx);
    } else {
	if (nx != sf_filesize(inp)) 
	    sf_error("Wrong model size");
    }

    model = sf_floatalloc(nx);
    
    if (adj) {
	sf_floatread(data,nd,inp);
	for (ix=0; ix < nx; ix++) {
	    model[ix] = 0.0f;
	}
    } else {
	sf_floatread(model,nx,inp);
	for (id=0; id < nd; id++) {
	    data[id] = 0.0f;
	}
    }

    sx=0;
    for (im=0; im < nm; im++) {
	for (id=0; id < nd; id++) {
	    ix = indx[im][id]+sx;
	    
	    if (adj) {
		model[ix] += data[id];
	    } else {
		data[id] += model[ix];
	    }
	}
	sx += size[im];
    }

    if (adj) {
	sf_floatwrite(model,nx,out);
    } else {	
	sf_floatwrite(data,nd,out);
    } 
    
    exit(0);
}
