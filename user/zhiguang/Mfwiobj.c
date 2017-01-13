/* Calculate the misfit fuction  in Full Waveform Inversion */

#include <rsf.h>

#include "fdprep.h"
#include "sparsesolver.h"
#include "optimization.h"
#include "waveoperator.h"

int main(int argc, char *argv[])
{
	int n1, n2, npml, pad1, pad2, ns, uts;
	float d1, d2, omega, misfit;
	sf_complex ***f, ***obs;
	sf_file in, out, source, receiver, record;
	char *order;
	float **v, **recloc;

	sf_init(argc, argv);
	
	in=sf_input("in");
	out=sf_output("out");
	source=sf_input("source");
	receiver=sf_input("receiver");
	record=sf_input("record");

	uts=1;

	if(!sf_getint("npml", &npml)) npml=20;
	if(!sf_getfloat("omega", &omega)) sf_error("Need input omega.");
	if(NULL==(order = sf_getstring("order"))) order="j";
	/* discretization scheme (default optimal 9-point) */

	fdprep_order(order);

	/* read input dimension */
	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input.");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input.");
	if (!sf_histfloat(in, "d1", &d1)) sf_error("No d1= in input.");
	if (!sf_histfloat(in, "d2", &d2)) sf_error("No d2= in input.");
	if (!sf_histint(record, "n3", &ns)) sf_error("No n3= in record.");

	/* set up output dimension */
	sf_putint(out, "n1", 1);
	sf_putint(out, "n2", 1);

	/* PML padding */
	pad1=n1+2*npml;
	pad2=n2+2*npml;

	v=sf_floatalloc2(n1, n2);
	sf_floatread(v[0], n1*n2, in);
	recloc=sf_floatalloc2(n1, n2);
	sf_floatread(recloc[0], n1*n2, receiver);

	f=sf_complexalloc3(n1, n2, ns);
	sf_complexread(f[0][0], n1*n2*ns, source);
	obs=sf_complexalloc3(n1, n2, ns);
	sf_complexread(obs[0][0], n1*n2*ns, record);

	misfit=forward_operator(uts, pad1, pad2, omega, n1, n2, d1, d2,
			npml, ns, f, obs, false, recloc, v);

	sf_floatwrite(&misfit, 1, out);

	exit(0);
}
