/* 1-D ENO interpolation */
#include <rsf.h>

int main(int argc, char* argv[])
{
    int   n1, nn1;
    float o1, oo1;
    float d1, dd1;
    int   i1, ii1;
    float f1;
    
    float *din;
    float *dou;

    float d,g;
    int order;
    
    sf_eno map1d;
    sf_file Fi, Fo, Fp;

    sf_init (argc, argv);
    Fi = sf_input ("in");
    Fo = sf_output("out");

    /*------------------------------------------------------------*/
    /* get input dimensions */
    if (!sf_histint  (Fi,"n1",&n1)) n1=1;
    if (!sf_histfloat(Fi,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(Fi,"d1",&d1)) d1=1.0;
    /*------------------------------------------------------------*/	
    if (NULL != sf_getstring("pattern")) {
	Fp = sf_input("pattern");
    } else {
	Fp = NULL;
    }
    /*------------------------------------------------------------*/
    /* get output dimensions */
    if(!sf_getint  ("n1",&nn1) && (NULL==Fp || !sf_histint  (Fp,"n1",&nn1))) nn1=n1;
    if(!sf_getfloat("d1",&dd1) && (NULL==Fp || !sf_histfloat(Fp,"d1",&dd1))) dd1=d1;
    if(!sf_getfloat("o1",&oo1) && (NULL==Fp || !sf_histfloat(Fp,"o1",&oo1))) oo1=o1;
    /*------------------------------------------------------------*/
    /* make output header */
    sf_putint  (Fo,"n1",nn1);
    sf_putfloat(Fo,"d1",dd1);
    sf_putfloat(Fo,"o1",oo1);
    /*------------------------------------------------------------*/        
    /* set-up interpolation */
    if (!sf_getint("order",&order)) order=3;
    order = SF_MIN( order,n1);
    map1d = sf_eno_init(order, n1);
    /*------------------------------------------------------------*/

    din = sf_floatalloc( n1);
    dou = sf_floatalloc(nn1);

    sf_floatread (din, n1,Fi);
    sf_eno_set(map1d,din);

    for(ii1=0; ii1<nn1; ii1++) { f1 = (oo1 + ii1 * dd1 - o1)/d1; i1 = f1; f1 -= i1;
	
	sf_eno_apply( map1d,
		      i1,
		      f1,
		      &d,&g,FUNC);
	dou[ii1] = d;
	
    }
    
    sf_floatwrite(dou,nn1,Fo);
    
    exit (0);
}

