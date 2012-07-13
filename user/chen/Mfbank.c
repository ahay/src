/* linear phase filter bank */

#include <rsf.h>
#include "dsp.h"

int main(int argc, char*argv[])
{
	sf_file in, out, flt;
	int n1, n2, n3, n[SF_MAX_DIM], axis, dim;
	int m1, m2, o1, d1, o2, d2;
	int i1, i2, i3, j2;
	float **wav, **fb, **c;
	char buf[8];

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");
	flt = sf_input("flt");

	dim = sf_filedims(in,n);
	if (!sf_histint(flt, "n1", &m1)) sf_error("No n1= in filter");
	if (!sf_histint(flt, "n2", &m2)) sf_error("No n2= in filter");
	if (!sf_histint(flt, "o1", &o1)) sf_error("No o1= in filter");
	if (!sf_histint(flt, "o2", &o2)) sf_error("No o2= in filter");
	if (!sf_histint(flt, "d1", &d1)) sf_error("No d1= in filter");
	if (!sf_histint(flt, "d2", &d2)) sf_error("No d2= in filter");

	if (!sf_getint("axis",&axis)) axis=1;
	/* which axis to apply the filter bank */

	if(axis>dim) sf_error("axis >= dim");

	n1 = 1; n2 = n[axis-1]; n3=1;
	for(i1=0; i1<axis-1; i1++ ) n1 *= n[i1];
	for(i1=axis; i1<dim; i1++ ) n3 *= n[i1];

	sf_shiftdim(in, out, axis);

	c = sf_floatalloc2(m1, m2);
	wav = sf_floatalloc2(n1, n2);
	fb  = sf_floatalloc2(n1, n2);

	sprintf(buf, "n%d", axis);
	sf_putint(out, buf, m2);
	sprintf(buf, "o%d", axis);
	sf_putfloat(out, buf, o2);
	sprintf(buf, "d%d", axis);
	sf_putfloat(out, buf, d2);

	sf_floatread(c[0], m1*m2, flt);

	for(i3=0; i3<n3; i3++)
	for(j2=0; j2<m2; j2++)
	{
		sf_floatread(wav[0], n1*n2, in);
		for(i1=0; i1<n1; i1++)
			firs(o1, o1+d1*(m1-1), c[j2]-o1, 
				wav[0]+i1, n1, n2, fb[0]+i1, n1);
		sf_floatwrite(fb[0], n1*n2, out);
	}

	free(wav[0]);
	free(wav);
	free(fb[0]);
	free(fb);
	free(c[0]);
	free(c);
	return 0;
}



