/* 2-D ENO interpolation */
#include <rsf.h>

int main(int argc, char* argv[])
{
    int   n1, n2, nn1, nn2;
    float o1, o2, oo1, oo2;
    float d1, d2, dd1, dd2;
    int   i1, i2, ii1, ii2;
    float f1, f2;
    
    float **din;
    float **dou;

    float d,g;
    int order;
    
    sf_eno2 map2d;
    sf_file Fi, Fo, Fp;

    sf_init (argc, argv);
    Fi = sf_input ("in");
    Fo = sf_output("out");

    /*------------------------------------------------------------*/
    /* get input dimensions */
    if (!sf_histint  (Fi,"n1",&n1)) n1=1;
    if (!sf_histfloat(Fi,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(Fi,"d1",&d1)) d1=1.0;

    if (!sf_histint  (Fi,"n2",&n2)) n2=1;
    if (!sf_histfloat(Fi,"o2",&o2)) o2=0.0;
    if (!sf_histfloat(Fi,"d2",&d2)) d2=1.0;
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

    if(!sf_getint  ("n2",&nn2) && (NULL==Fp || !sf_histint  (Fp,"n2",&nn2))) nn2=n2;
    if(!sf_getfloat("d2",&dd2) && (NULL==Fp || !sf_histfloat(Fp,"d2",&dd2))) dd2=d2;
    if(!sf_getfloat("o2",&oo2) && (NULL==Fp || !sf_histfloat(Fp,"o2",&oo2))) oo2=o2;
    /*------------------------------------------------------------*/
    /* make output header */
    sf_putint  (Fo,"n1",nn1);
    sf_putfloat(Fo,"d1",dd1);
    sf_putfloat(Fo,"o1",oo1);

    sf_putint  (Fo,"n2",nn2);
    sf_putfloat(Fo,"d2",dd2);
    sf_putfloat(Fo,"o2",oo2);
    /*------------------------------------------------------------*/        
    /* set-up interpolation */
    if (!sf_getint("order",&order)) order=3;
    order = SF_MIN(SF_MIN( order,n1),n2);
    map2d = sf_eno2_init(order, n1,n2);
    sf_warning("A");
    /*------------------------------------------------------------*/

    din = sf_floatalloc2( n1, n2);
    dou = sf_floatalloc2(nn1,nn2);

    sf_floatread (din[0], n1* n2,Fi);
    sf_eno2_set(map2d,din);

    for    (ii2=0; ii2<nn2; ii2++) { f2 = (oo2 + ii2 * dd2 - o2)/d2; i2 = f2; f2 -= i2;
	for(ii1=0; ii1<nn1; ii1++) { f1 = (oo1 + ii1 * dd1 - o1)/d1; i1 = f1; f1 -= i1;
	    
	    sf_eno2_apply( map2d,
			   i1,i2,
			   f1,f2,
			   &d,&g,FUNC);
	    dou[ii2][ii1] = d;
		
	}
    }
    
    sf_floatwrite(dou[0],nn1*nn2,Fo);
    
    exit (0);
}

