/* Accumulate */
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define CICLOOP2D(a)				\
    for    (i2=0;i2<sf_n(a2);i2++){		\
	for(i1=0;i1<sf_n(a1);i1++){		\
            {a}                                 \
        }                                       \
    }

#define CICLOOP3D(a)                                  \
    for        (i3=0;i3<sf_n(a3);i3++){		      \
	for    (i2=0;i2<sf_n(a2);i2++){		      \
            for(i1=0;i1<sf_n(a1);i1++){		      \
                {a}				      \
            }					      \
        }					      \
    }

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    bool verb;     /* verbosity flag */
    int  axis;

    sf_file  Fin, Fou; /* I/O files */

    float   **in2d=NULL, **ou2d=NULL;
    float  ***in3d=NULL,***ou3d=NULL;

    sf_axis a1,a2,a3,a4=NULL,aa; /* cube axes */
    int     i1,i2,i3,i4;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint( "axis",&axis)) axis=3; /* accumulate axis */
    
    /*------------------------------------------------------------*/
    Fin = sf_input ("in" );
    Fou = sf_output("out");

    a1=sf_iaxa(Fin,1); if(verb) sf_raxa(a1); 
    a2=sf_iaxa(Fin,2); if(verb) sf_raxa(a2); 
    a3=sf_iaxa(Fin,3); if(verb) sf_raxa(a3); 
    if(axis>3) 
	a4=sf_iaxa(Fin,4); if(verb) sf_raxa(a4); 
    aa=sf_maxa(1,0,1); sf_setlabel(aa,"");  sf_setunit (aa,""); 
    
    /*------------------------------------------------------------*/
    if(axis>3) {

    	in3d = sf_floatalloc3(sf_n(a1),sf_n(a2),sf_n(a3)); 
	ou3d = sf_floatalloc3(sf_n(a1),sf_n(a2),sf_n(a3));
	CICLOOP3D( ou3d[i3][i2][i1]=0; );

	for(i4=0; i4<sf_n(a4); i4++) {
	    sf_floatread(in3d[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fin);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(i1,i2,i3)							\
    shared (a1,a2,a3,in3d,ou3d)
#endif
	    CICLOOP3D( ou3d[i3][i2][i1] += in3d[i3][i2][i1]; );
	    sf_floatwrite(ou3d[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fou); 
	}
	
	free(**in3d); free(*in3d); free(in3d);
	free(**ou3d); free(*ou3d); free(ou3d);

    } else {

	in2d = sf_floatalloc2(sf_n(a1),sf_n(a2)); 
	ou2d = sf_floatalloc2(sf_n(a1),sf_n(a2));
	CICLOOP2D( ou2d[i2][i1]=0.; );

	for(i3=0; i3<sf_n(a3); i3++) {
	    sf_floatread(in2d[0],sf_n(a1)*sf_n(a2),Fin);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(i1,i2)							\
    shared (a1,a2,in2d,ou2d)
#endif
	    CICLOOP2D( ou2d[i2][i1] += in2d[i2][i1]; );
	    sf_floatwrite(ou2d[0],sf_n(a1)*sf_n(a2),Fou); 
	}
	
	free(*in2d); free(in2d);
	free(*ou2d); free(ou2d);
    }

    exit (0);
}
