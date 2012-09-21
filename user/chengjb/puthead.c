/*************************************************************************
* *  put head for rsf datasets
* *
* *    Copyright: Tongji University (Jiubing Cheng)
* *    2012.3.2
* *************************************************************************/
#include <rsf.h>

void puthead1(sf_file Fo, int n1, int n2, float d1, float o1, float d2, float o2)
/*<  put head for (kx,kz) domain float-type data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","kz");
        sf_putstring(Fo,"label2","kx");
        sf_putstring(Fo,"unit1","2*pi/m");
        sf_putstring(Fo,"unit2","2*pi/m");
}

void puthead2(sf_file Fo, int n1, int n2, float d1, float o1, float d2, float o2)
/*<  put head for (x,z) domain float-type data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"unit1","km");
        sf_putstring(Fo,"unit2","km");
}

void puthead3(sf_file Fo, int n1, int n2, int n3, float d1, float d2, float d3, float o1, float o2)
/*<  put head for (x,z,t) domain float-type 3D data sets>*/
{
        /* Read/Write axes */
        sf_putint(Fo,"n1",n1);
        sf_putint(Fo,"n2",n2);
        sf_putint(Fo,"n3",n3);
        sf_putfloat(Fo,"d1",d1);
        sf_putfloat(Fo,"d2",d2);
        sf_putfloat(Fo,"d3",d3);
        sf_putfloat(Fo,"o1",o1);
        sf_putfloat(Fo,"o2",o2);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"label3","t");
        sf_putstring(Fo,"unit1","km");
        sf_putstring(Fo,"unit2","km");
        sf_putstring(Fo,"unit3","second");
}
