/* Ormsby wavelet (sfscale axis=1 can be used to normalize the function) */

#include <rsf.h>
#include <math.h>

int main(int argc, char *argv[])
{
	int nt, it;
	float f1, f2, f3, f4;
	float dt, t0, shift;
	float *f, t, tmp1, tmp2, tmp3, tmp4;
	sf_file Fout;

	sf_init(argc, argv);

	Fout=sf_output("out");
	if(!sf_getint("nt", &nt)) sf_error("Need to input nt=");
	if(!sf_getfloat("dt", &dt)) dt=0.002;
	if(!sf_getfloat("t0", &t0)) t0=0.;

	if(!sf_getfloat("f1", &f1)) f1=5.0;
	/* low-cut frequency */
	if(!sf_getfloat("f2", &f2)) f2=10.0;
	/* low-pass frequency */
	if(!sf_getfloat("f3", &f3)) f3=55.0;
	/* high-pass frequency */
	if(!sf_getfloat("f4", &f4)) f4=60.0;
	/* high-cut frequency */
	if(!sf_getfloat("shift", &shift)) shift=0.1;

	f=sf_floatalloc(nt);
	sf_putint(Fout, "n1", nt);
	sf_putfloat(Fout, "d1", dt);
	sf_putfloat(Fout, "o1", t0);

	for(it=0; it<nt; it++){
		t=it*dt+t0-shift;
		tmp1=f1*SF_PI;
		tmp2=f2*SF_PI;
		tmp3=f3*SF_PI;
		tmp4=f4*SF_PI;
		if(t==0.){
			f[it]=tmp4+tmp3-tmp2-tmp1;
		}else{
			tmp4=(sin(tmp4*t)*sin(tmp4*t)-sin(tmp3*t)*sin(tmp3*t))/(tmp4-tmp3);
			tmp2=(sin(tmp2*t)*sin(tmp2*t)-sin(tmp1*t)*sin(tmp1*t))/(tmp2-tmp1);
			f[it]=(tmp4-tmp2)/t/t;
		}
	}
	sf_floatwrite(f, nt, Fout);
	sf_fileclose(Fout);
	exit(0);
}
