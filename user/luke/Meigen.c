

#include <rsf.h>


// make a vector from data file
void makevec(float *vec, float** data, int i, int nx){
	int j;
	for( j=0; j < nx; j++){
		vec[j] = data[j][i];				
	}	
	return;	
}

void dotproduct(float *dot, float *vec1, float *vec2, int nx){
	int i;
	*dot = 0.0; //initialize
	for (i=0; i<nx ; i++){
		*dot += vec1[i]*vec2[i];
	}
	return;
}

int main(int argc, char *argv[])
{	//declare variables
	int i,ei;
	int nt, nx;

	float dt, dx, t0, x0;


	float **matrix=NULL;
	float *eigen=NULL;
\

	sf_file in=NULL, out=NULL;
	sf_init (argc,argv);
	in = sf_input("in");
	out = sf_output("out");

	//read data
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");


	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");


	if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o1= in input");

	// allocate arrays
	matrix = sf_floatalloc2(nt,nx);
	eigen = sf_floatalloc(nt);


	// setup output
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);

	sf_putint(out,"n2",1);
	sf_putfloat(out,"d2",1);
	sf_putfloat(out,"o2",0);


	// read data
	sf_floatread (matrix[0],nt*nx,in);

	ei = 0;
//ACTUALLY CALL PROGRAM
	int N = nx;
	int LDA = N;
	float *WR;
	float *WI;
	float *VL;
	int LDVL = N;
	float **VR;
	int LDVR = N;
	float *WORK;
	int LWORK = 3*N;
	int INFO;
	// dsytrd_(matrix,5,matrix,matrix);
	//sgeev_(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
	sgeev_('N','V',N,matrix,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO);
	
	// write data
	sf_floatwrite(matrix,nt,out);
	exit(0);
}
