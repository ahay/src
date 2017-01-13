/* Phase encoding */

#include <rsf.h>
#include <time.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	bool encoding;
	int i1, i2, i3, n1, n2;
	int ns, nsim, nsource, belong;
	float *random;
	sf_complex phase;
	sf_complex ***newarray, **oldarray;
	sf_file in, out, oldrec, newrec;

	sf_init(argc, argv);

	/* I/O */
	in=sf_input("in");
	out=sf_output("out");
	oldrec=sf_input("oldrec");
	newrec=sf_output("newrec");

	/* read arguments from file */
	if(!sf_histint(in, "n1", &n1)) sf_error("No n1 in input.");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2 in input.");
	if(!sf_histint(in, "n3", &ns)) sf_error("No n3 in input.");

	/* read parameters from command-line */
	if(!sf_getbool("encoding", &encoding)) encoding=true;
	if(!sf_getint("nsim", &nsim)) nsim=ns;
	if(!sf_getint("nsource", &nsource)) nsource=1;

	/* check */
	if((ns-1)/nsource+1 != nsim) sf_error("ns/nsource should be equal to nsim!");

	/* set up the dimension of the output files */
	sf_putint(out, "n3", nsim);
	sf_putint(newrec, "n3", nsim);

	/* storage allocation */
	newarray=sf_complexalloc3(n1, n2, nsim);
	oldarray=sf_complexalloc2(n1, n2);
	random=sf_floatalloc(ns);

	/* generate random number */
	srand(time(NULL));
	for(i1=0; i1<ns; i1++)
		random[i1]=((double)(2.*3.1415926)/RAND_MAX)*rand();

	/* Initialization */
	for(i3=0; i3<nsim; i3++)
		for(i2=0; i2<n2; i2++)
			for(i1=0; i1<n1; i1++)
				newarray[i3][i2][i1]=sf_cmplx(0.,0.);

	/* Big Loop for file 1*/
	for(i3=0; i3<ns; i3++){
		if(encoding)
			phase=sf_cmplx(cosf(random[i3]), sinf(random[i3]));
		else
			phase=sf_cmplx(1.,0.);
		belong=i3%nsim;
		
		sf_complexread(oldarray[0], n1*n2, in);

		for(i2=0; i2<n2; i2++){
			for(i1=0; i1<n1; i1++){
#ifdef SF_HAS_COMPLEX_H
				oldarray[i2][i1] = oldarray[i2][i1]*phase;
#else
				oldarray[i2][i1] = sf_cmul(oldarray[i2][i1], phase);
#endif

#ifdef SF_HAS_COMPLEX_H
				newarray[belong][i2][i1] += oldarray[i2][i1];
#else
				newarray[belong][i2][i1] = sf_cadd(newarray[belong][i2][i1], oldarray[i2][i1]);
#endif
			}
		}
	}

	sf_complexwrite(newarray[0][0], n1*n2*nsim, out);

	/* Initialization */
	for(i3=0; i3<nsim; i3++)
		for(i2=0; i2<n2; i2++)
			for(i1=0; i1<n1; i1++)
				newarray[i3][i2][i1]=sf_cmplx(0.,0.);

	/* Big Loop for file 2*/
	for(i3=0; i3<ns; i3++){
		if(encoding)
			phase=sf_cmplx(cosf(random[i3]), sinf(random[i3]));
		else
			phase=sf_cmplx(1.,0.);
		belong=i3%nsim;
		
		sf_complexread(oldarray[0], n1*n2, oldrec);

		for(i2=0; i2<n2; i2++){
			for(i1=0; i1<n1; i1++){
#ifdef SF_HAS_COMPLEX_H
				oldarray[i2][i1] = oldarray[i2][i1]*phase;
#else
				oldarray[i2][i1] = sf_cmul(oldarray[i2][i1], phase);
#endif

#ifdef SF_HAS_COMPLEX_H
				newarray[belong][i2][i1] += oldarray[i2][i1];
#else
				newarray[belong][i2][i1] = sf_cadd(newarray[belong][i2][i1], oldarray[i2][i1]);
#endif
			}
		}
	}

	sf_complexwrite(newarray[0][0], n1*n2*nsim, newrec);

	exit(0);
}
