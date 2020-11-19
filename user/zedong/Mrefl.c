/* Generate the reflector which will be used in PERM. The input is: 
 *velocity:nh,nx,nz
 * for the slice h=0,refl(h=0,x,z)=laplace(velocity(h=0,x,z))
 * refl(h!=0,x,z)=0
 * */
/*uploaded by zedong wu*/
#include <rsf.h>
#define LONG_INT long long 
#define LEN_FILE 3000
#define _FILE_OFFSET_BITS 64
int main(int argc, char* argv[])
{
    sf_init(argc,argv);
    sf_file Fin,Fout;
    long long int n1,n2,n3;
    float d1,o1;
    Fin=sf_input ("in" );
    sf_axis a1,a2,a3;
    a1 = sf_iaxa(Fin,1); n1 = sf_n(a1); d1=sf_d(a1);o1=sf_o(a1); 
    a2=sf_iaxa(Fin,2);n2=sf_n(a2);
    a3=sf_iaxa(Fin,3);n3=sf_n(a3);
    fprintf(stderr,"%lld %lld %lld %f %f\n",n1,n2,n3,d1,o1);
    long long int m=((long long int)n1)*n2*n3;
    float *vels=(float *)malloc(sizeof(float)*m);
    char* filename=sf_histstring(Fin,"in");
    FILE *fp=fopen(filename,"rb");
    if (!fread(vels,sizeof(float),m,fp))
        abort();
    fclose(fp);
    int i2,i3,i1;
    //fp=fopen("a.dat","wb");
    float *refl=(float *)calloc(m,sizeof(float));
    int icenter=0;
    float offset=1e10;
    for(i1=0;i1<n1;i1++)
    {
        if(fabs(o1+d1*i1)<offset)
        {
            offset=fabs(o1+d1*i1);
            icenter=i1;
         
	}
    }
    fprintf(stderr,"minoffset=%f %d\n",offset,icenter);
    float max=0;
    for(i3=1;i3<(n3-1);i3++)
    for(i2=1;i2<(n2-1);i2++)
    for(i1=0;i1<n1;i1++)
    {
        long long int ij=i3*n2*n1+i2*n1+i1;
        float a=vels[ij-n2*n1]+vels[ij+n2*n1]+vels[ij-n1]+vels[ij+n1]-4*vels[ij];
        if(i1==icenter)
        {
	    refl[ij]=a;
	    if(a>max)
	    max=a;
	}
    }
    fprintf(stderr,"maxvalue=%22.16f\n",max);
   // for(i3=0;i3<n3;i3++)
   // for(i2=0;i2<n2;i2++)
   // {
    //    long long int ij=i3*n2*n1+i2*n1+icenter;
    //    fwrite(refl+ij,sizeof(float),1,fp);
   // }
    //fclose(fp);
    Fout = sf_output("out"); 
    sf_oaxa(Fout,a1,1);
    sf_oaxa(Fout,a2,2);
    sf_oaxa(Fout,a3,3);
    sf_floatwrite(refl,m,Fout);
    sf_close();
    return 0;
}
